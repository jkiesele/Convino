
The Convino program allows combining sets of differential or inclusive measurements including but not limited to measurements obtained with simultaneously fitting a set of nuisance parameters, representing sources of correlated systematic uncertainties. In these fits these nuisance parameters are typically constrained by the prior knowledge on the systematic uncertainties as well as by the data. As a result of this procedure all fitted parameters are correlated among each other. The best approach for a combination of these measurements would be the maximization of a combined likelihood. To define this likelihood, the full fit model used for each measurement and the original data are required. However, only in rare cases this information is publicly available. In absence of this information most commonly used combinations methods are not able to account for such correlations between uncertainties. The method used by the Convino program provides a solution for this problem. It relies on the public result and its covariance or Hessian, only, and was validated against the combined-likelihood approach. This dedicated software package provides a text-based user interface alongside a C++ interface. The latter integrates ROOT histogram and graph classes for simple combination of binned measurements such as differential cross sections.

For a description of the method and the program, see [1]. Please also cite this document if this tool is used.
New features with respect to the ones described in the published paper are described in the file documentation.pdf, wich is provided in this repository.

The source code can be compiled with gcc, supporting C++11 standard. It depends on ROOT6 [https://root.cern.ch/].

It is also installed on lxplus6 at CERN in
``/afs/cern.ch/user/j/jkiesele/public/Convino/latest``
For compiling and running on cern lxplus, please source the lxplus_env.sh script.


[1] J. Kieseler, "A method and tool for combining differential or inclusive measurements obtained with simultaneously constrained uncertainties", https://arxiv.org/abs/1706.01681, EPJC (2017) 77
