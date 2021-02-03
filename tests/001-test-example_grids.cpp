//
// example grids
//

#include "AMReX.H"
#include "AMReX_BaseFab.H"
#include "AMReX_BoxArray.H"
#include "AMReX_FArrayBox.H"
#include "AMReX_MultiFab.H"

#include "AMReX_DistributionMapping.H"
#include "AMReX_Geometry.H"
#include "AMReX_ParmParse.H"
#include "AMReX_Print.H"

void print_int_vec(const amrex::IntVect &iv) {
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        amrex::Print() << "iv[" << i << "] = " << iv[i] << '\n';
    }
}

void example_int_vec() {
    // play with int vec
    amrex::IntVect iv(AMREX_D_DECL(19, 0, 5));
    print_int_vec(iv);

    amrex::IntVect nv(AMREX_D_DECL(127, 127, 127));
    amrex::IntVect coarsening_ratio(AMREX_D_DECL(2, 2, 2));
    nv.coarsen(2);
    print_int_vec(nv);
    nv.coarsen(coarsening_ratio);
    print_int_vec(nv);

    amrex::IntVect v = amrex::coarsen(nv, coarsening_ratio);
    print_int_vec(v);
}

void example_index_type() {
    amrex::IndexType xface(amrex::IntVect{AMREX_D_DECL(1, 0, 0)});

    amrex::Print() << "Is cell centred: " << xface.cellCentered() << '\n';
}

void example_box() {
    amrex::IntVect lower{AMREX_D_DECL(64, 64, 64)};
    amrex::IntVect upper{AMREX_D_DECL(127, 127, 127)};

    // node based in x,y and z
    amrex::IndexType type({AMREX_D_DECL(1, 1, 1)});

    // nodal box
    amrex::Box cb(lower, upper);
    amrex::Box nd(lower, upper + 1, type);

    amrex::Print() << "Cell box cb: " << cb << '\n';
    amrex::Print() << "Nodal box nd: " << nd << '\n';

    // A new Box with type (node, node, node)
    amrex::Box nd1 = amrex::surroundingNodes(cb);
    amrex::Print() << "A new box nd1: " << nd1 << '\n';

    amrex::Box cb1 = amrex::enclosedCells(nd);
    amrex::Print() << "A new box cb1: " << cb1 << '\n';
    if (cb == cb1) {
        amrex::Print() << "Box cb and cb1 are equal\n";
    }

    amrex::IndexType new_type({AMREX_D_DECL(0, 1, 1)});
    amrex::Box b3 = amrex::convert(cb, new_type);

    amrex::Print() << "A new box b3: " << b3 << '\n';
    amrex::Print() << "Surrounding nodes of b3: " << b3.surroundingNodes() << '\n';
    amrex::Print() << "Enclosed cells of b3: " << b3.enclosedCells() << '\n';

    // access internal data of a box:
    const auto &small_end = cb.smallEnd();
    const auto &big_end = cb.bigEnd();
    amrex::Print() << "Small/big ends of box cb: " << small_end << ", " << big_end << '\n';

    // refine and coarse a box
    // the box covers always the same space!
    cb.refine(2);
    amrex::Print() << "cb fine: " << cb << '\n';
    cb.coarsen(2);
    amrex::Print() << "cb coarse: " << cb << '\n';

    // loop over a box:
    const auto &lo = amrex::lbound(cb);
    amrex::Print() << "lower bound of cb: " << lo << '\n';
    const auto &hi = amrex::ubound(cb);
    amrex::Print() << "upper bound of cb: " << hi << '\n';
    const auto &box_length = amrex::length(cb);
    amrex::Print() << "box length of cb: " << box_length << '\n';

    for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                //amrex::Print() << "(" << i << ',' << j << ',' << k << ")\n";
            }
        }
    }

    // real box and geometry
    // coord type: 0 cartesian, 1 cylindrical and 2 spherical
    // bc type: 0 non-periodic, 1 periodic

    constexpr int n_node = 64;
    amrex::Box domain(amrex::IntVect{AMREX_D_DECL(0, 0, 0)},
                      amrex::IntVect{AMREX_D_DECL(n_node, n_node, n_node)},
                      type);


    amrex::RealBox real_box({AMREX_D_DECL(-1, -1, -1)}, {AMREX_D_DECL(1, 1, 1)});
    constexpr int coord_sys = 0;

    std::array<int, AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(1, 1, 1)};

    amrex::Geometry geom(domain, real_box, coord_sys, is_periodic);
    amrex::Print() << "Geometry: " << geom << '\n';
}

void example_box_array() {
    // make a new box (cell,cell,cell)
    amrex::Box domain(amrex::IntVect{AMREX_D_DECL(0, 0, 0)}, amrex::IntVect{AMREX_D_DECL(127, 127, 127)});

    amrex::BoxArray ba(domain);
    amrex::Print() << "BoxArray size is " << ba.size() << "\n";// 1
    ba.maxSize(64);                                            // Chop into boxes of 64^spacedim cells
    amrex::Print() << ba;

    amrex::DistributionMapping dm(ba);
    amrex::Print() << "Distribution mapping: " << dm << '\n';
    //
    // amrex::Vector<int> pmap {...};
    // The user fills the pmap array with the values specifying owner processes
    // dm.define(pmap);  // Build DistributionMapping given an array of process IDs.
}

void example_base_fab() {
    amrex::Box bx(amrex::IntVect{AMREX_D_DECL(-4, 8, 32)},
                  amrex::IntVect{AMREX_D_DECL(32, 64, 48)});
    constexpr int num_comps = 4;

    // data field on box with num_comps
    amrex::BaseFab<amrex::Real> fab(bx, num_comps);

    // get the box back where we have data field fab;
    const amrex::Box &box_of_fab = fab.box();

    amrex::Print() << "Comps of fab: " << fab.nComp() << '\n';

    // get a pointer to thw array data
    // to nth component
    auto data_ptr = fab.dataPtr(0);
    amrex::Print() << "pointer to nth comp data: " << *data_ptr << '\n';

    auto min_val = fab.min(2);
    amrex::Print() << "min value of comp 2: " << min_val << '\n';

    amrex::FArrayBox fab1(bx, num_comps);
    amrex::FArrayBox fab2(bx, num_comps);

    fab1.setVal(1.0);
    fab1.mult(10.0, 0);
    fab2.setVal(2.0);
    amrex::Real a = 3.0;
    fab2.saxpy(a, fab1);// for all comps fab2 <- a*fab1 + fab2
    amrex::Print() << "fab2: " << fab2 << '\n';

    // to define if function should be called on host (cpu) or device (gpu)
    // if amrex is built without gpu backend, device becomes host (cpu)
    fab1.setVal<amrex::RunOn::Host>(3.0);
    fab1.mult<amrex::RunOn::Device>(10.0, 3);

    // for performing data on array use:
    const amrex::Array4<amrex::Real> &fab1_arr = fab1.array();
    const amrex::Array4<const amrex::Real> fab2_arr = fab2.const_array();
}

void example_multi_fab() {
    amrex::Box bx(amrex::IntVect{AMREX_D_DECL(0, 0, 0)},
                  amrex::IntVect{AMREX_D_DECL(16, 16, 16)});

    amrex::Print() << "Box: " << bx << '\n';

    amrex::BoxArray ba(bx);

    amrex::DistributionMapping dm(ba);

    int ncomp = 4; // number of components
    int ngrow = 1; // number of ghost layers

    // cell based
    amrex::MultiFab mf(ba, dm, ncomp, ngrow);

    
    // node based multi fab on same box array
    amrex::MultiFab mf_node(amrex::convert(ba, amrex::IntVect{AMREX_D_DECL(1,1,1)}), dm, ncomp, ngrow);
    amrex::Print() << "Multi fab: " << mf.size() << '\n';
}

int main(int argc, char *argv[]) {
    amrex::Initialize(argc, argv);
    amrex::Print() << "Hello world from AMReX version "
                   << amrex::Version() << "\n";

    amrex::Print() << "Dimension: " << AMREX_SPACEDIM << '\n';

    //example_int_vec();

    //example_index_type();

    //example_box();

    //example_box_array();

    //example_base_fab();

    example_multi_fab();

    amrex::Finalize();
}