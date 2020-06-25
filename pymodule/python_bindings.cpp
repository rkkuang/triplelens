
// #include <pybind11/pybind11.h>
// #include <pybind11/stl.h>
#include </Library/Python/3.7m/include/pybind11/pybind11.h>
#include </Library/Python/3.7m/include/pybind11/stl.h>
#include "../src/VBBinaryLensingLibrary.h"
#include "../src/TripleLensingLibrary.h"
#include <string>

namespace py = pybind11;


TripleLensing TRIL;
// // Declaration of an instance to VBBinaryLensing class. 
// VBBinaryLensing VBBL;

PYBIND11_MODULE(TripleLensing, m) {
    py::options options;
    options.disable_function_signatures();
    py::class_<TripleLensing>(m, "TripleLensing")
        .def(py::init())
        // Settings
        // .def("LoadESPLTable", &VBBinaryLensing::LoadESPLTable,
            // """Loads a pre calculated binary table for extended source calculation.""")
        .def_readwrite("secnum", &TripleLensing::secnum,
                "divide source boundary into how many parts.")
        .def_readwrite("basenum", &TripleLensing::basenum,
                "how many points in each part.")
        .def_readwrite("nphi", &TripleLensing::nphi,
                "how many points to start.")
        .def_readwrite("quaderr_Tol", &TripleLensing::quaderr_Tol,
                "Tolerance of quadrupole test.")
        // Maginfication calculations
        // .def("tripleFS2python", &TripleLensing::tripleFS2python,
        //     py::return_value_policy::reference,
        //     R"mydelimiter(
        //     Extended Source Point Lens magnification calculation.

        //     Magnification of a uniform brightness-source by a triple lens system.
        //     Parameters
        //     ----------
        //     u : float 
        //         Distance of source from the center of the lens.
        //     rho : float 
        //         Source radius in units of the Einstein radius of the lens.

        //     Returns
        //     -------
        //     float
        //         Magnification.
        //     )mydelimiter")

        .def("tripleFS2python", 
            [](TripleLensing &self, std::vector<double> mlens, std::vector<double> zlens, double xsCenter, double ysCenter, double rs){
                double mu;
                mu = self.tripleFS2python(mlens.data(), zlens.data(), xsCenter, ysCenter, rs);
                return mu; 
            },
            R"mydelimiter(
            outputCriticalTriple_list
            )mydelimiter")
        .def("solv_lens_equation", 
            [](TripleLensing &self, std::vector<double> mlens, std::vector<double> zlens,
                                      double xs, double ys, int nlens )
            {
                std::vector<double> zrxy( (nlens*nlens+1)*2 );
                self.solv_lens_equation(zrxy.data(), mlens.data(), zlens.data(),xs,ys,nlens
                        );
                return zrxy;
            },
            R"mydelimiter(
            solv_lens_equation
            )mydelimiter")
        .def("outputCriticalTriple_list", 
            [](TripleLensing &self, std::vector<double> mlens, std::vector<double> zlens, int nlens, int NPS )
            {
                std::vector<double> allxys(NPS*40);
                // std::vector<double> allxys(NPS*4+2);
                self.outputCriticalTriple_list(allxys.data(), mlens.data(), zlens.data(),nlens,NPS);
                return allxys;
            },
            R"mydelimiter(
            outputCriticalTriple_list
            )mydelimiter");

//     py::class_<_theta>(m, "_theta")
//         .def(py::init<double>()); //constructor 

//     py::class_<_point>(m, "_point")
//         .def(py::init<double, double, _theta *>()) 
//         .def_readwrite("next", &_point::next)
//         .def_readwrite("prev", &_point::prev)
//         .def_readonly("x1", &_point::x1)
//         .def_readonly("x2", &_point::x2);

//     py::class_<_curve>(m, "_curve")
//         .def(py::init<_point*>()) //constructor 1
//         .def(py::init()) //constructor 2
//         .def_readwrite("first", &_curve::first)
//         .def_readwrite("last", &_curve::last)
//         .def_readwrite("next", &_curve::next)
//         .def_readwrite("prev", &_curve::prev);

//     py::class_<_sols>(m, "_sols")
//         .def(py::init()) //constructor
//         .def_readwrite("first", &_sols::first)
//         .def_readwrite("last", &_sols::last);

}