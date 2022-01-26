#include <ctime>
#include <fstream>
#include <iostream>
#include <vector>
#include <streambuf>
#include <string>
#include <thread>

#include "math/tiny/tiny_double_utils.h"

#include "math/tiny/tiny_dual.h"
#include "math/tiny/tiny_dual_double_utils.h"

#include "math/tiny/tiny_algebra.hpp"

#include "utils/file_utils.hpp"
#include "urdf/urdf_parser.hpp"
#include "urdf/urdf_to_multi_body.hpp"
#include "dynamics/forward_dynamics.hpp"
#include "dynamics/integrator.hpp"
//#include "math/tiny/tiny_float_utils.h"

#include "urdf/urdf_cache.hpp"

using namespace TINY;
using namespace tds;

// typedef double TinyDualScalar;
typedef double InnerScalar;
typedef ::TINY::TinyDualDouble MyScalar;
typedef ::TINY::TinyDualDoubleUtils MyTinyConstants;

typedef TinyAlgebra<MyScalar, MyTinyConstants> MyAlgebra;
typedef tds::MultiBodyContactPoint<MyAlgebra> MultiBodyContactPoint;

// MyAlgebra::Vector3 gravity( MyAlgebra::zero(),
//                             MyAlgebra::zero(),
//                             MyAlgebra::fraction(-980, 100));

 auto ts = MyAlgebra::fraction(1, 100);

int main(){
    std::string anymal_path = "/home/sgiles/anymal_ws/src/anymal_c_simple_description/urdf/anymal2.urdf";
    std::string plane_path = "/home/sgiles/src/tiny-differentiable-simulator/data/plane_implicit.urdf";
    UrdfParser<MyAlgebra> parser;
    UrdfStructures urdf_data = parser.load_urdf(anymal_path);
    UrdfStructures plane_urdf_data = parser.load_urdf(plane_path);
    
    MultiBody<MyAlgebra> plane_mb(false);
    MultiBody<MyAlgebra> mb(true);

    World<MyAlgebra> world;

    UrdfToMultiBody<MyAlgebra>::convert_to_multi_body(plane_urdf_data, world, plane_mb, 0);
    plane_mb.initialize();

    UrdfToMultiBody<MyAlgebra>::convert_to_multi_body(urdf_data, world, mb, 0);
    mb.initialize();

    std::vector<MultiBody<MyAlgebra> *> bodies;
    bodies.emplace_back(&plane_mb);
    bodies.emplace_back(&mb);
    auto dispatcher = world.get_collision_dispatcher();
    MultiBodyConstraintSolver<MyAlgebra> mb_solver;
    mb_solver.pgs_iterations_ = 50;

    // std::cout << "Initialized" << std::endl;
    // std::cout << "Enter to run.";
    // std::cin.get();

    mb.q_[6].set_real(0.65);
    mb.tau_[0].set_dual(1.);

    auto qdd_out = mb.qd_;


    auto contacts = world.compute_contacts_multi_body(bodies, &dispatcher);
    auto collisions = mb_solver.resolve_collision2(contacts[0], ts);

    int N = 1000;

    std::clock_t c_start = std::clock();
    for (int i=0; i < N; i++){
        for (int j=0; j < 18; j++){
            mb.qd_[j].set_real(0.0);
            mb.qd_[j].set_dual(0.0);
        }
        mb.qd_[i % 18].set_dual(1.0);


        auto qd_out = mb.qd_;

        forward_kinematics(mb, mb.q(), mb.qd());
        // Assign input to mb.tau_


        forward_dynamics(mb, world.get_gravity());
        integrate_euler_qdd(mb, ts);

        contacts = world.compute_contacts_multi_body(bodies, &dispatcher);
        collisions = mb_solver.resolve_collision2(contacts[0], ts);
        // qdd_out = [(mb.qd_ - qd_out) / ts for i in range(mb.qd.size())]
        qdd_out = (mb.qd_ - qd_out) / ts;
                
    }
    std::clock_t c_end = std::clock();
    std::cout << qdd_out[0].real() << '\n';

    // for (int i=0; i<18; i++){
    //     auto x = qd_out[i];
    //     std::cout << x.real() << '\t' << x.dual() << std::endl;
    // }    
    // for (int i=0; i<18; i++){
    //     auto x = qdd_out[i];
    //     std::cout << x.real() << '\t' << x.dual() << std::endl;
    // }

    // integrate_euler(mb, ts);


    double time_elapsed_ms = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC / N;
    std::cout << "collisions: " << collisions.size() << "\n";
    std::cout << "CPU time used: " << time_elapsed_ms << " ms\n";
    std::cout << "CPU time used: " << c_end-c_start << " clocks\n";
    return 0;
}