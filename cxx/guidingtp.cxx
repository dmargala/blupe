// Created 18-Jul-2013 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>
// convert guiding offsets to tp correction 
// g++ -lboost_program_options guidingtp.cxx -o guidingtp

#include "boost/program_options.hpp"
#include "boost/math/special_functions/gamma.hpp"

#include <iostream>
#include <fstream>
#include <cmath>

namespace po = boost::program_options;


const double PLATESCALE = 217.7358/3600; /* mm/arcsec */

// Calculates the fraction of light from an object with the specified Gaussian
// fwhm that enters a fiber of the specified diameter when the object's centroid
// is offset from the fiber center by the specified amount. All inputs should
// be specified in the same units. The default diameter is the BOSS fiber size
// in arcsecs, so requires that fwhm and offset also be specified in arcsecs.
// Accuracy gets worse with increasing diameter/fwhm, but should be at least 0.5%
// for fwhm > diameter/4 (or fwhm > 0.5 arcsec for BOSS fibers).
double fiberFraction(double fwhm, double offset, double diameter=2.0/*arcsec*/) {
    if(fwhm <= 0) {
        std::cerr << "fiberFraction: invalid fwhm <= 0." << std::endl;
    }
    if(offset < 0) {
        std::cerr << "fiberFraction: invalid offset < 0." << std::endl;
    }
    if(diameter <= 0) {
        std::cerr << "fiberFraction: invalid diameter <= 0." << std::endl;
    }
    if(diameter > 4*fwhm) {
        std::cerr << "fiberFraction: diameter > 4*fwhm not implemented yet." << std::endl;
    }
    // Convert from FWHM to sigma.
    double sigma(fwhm/2.3548200450309493820); // constant is 2*sqrt(2*log(2))
    // Calculate dimensionless ratios
    double t(offset/sigma), ss(diameter/(2*sigma));
    double tSqby2(t*t/2),ssSqby2(ss*ss/2);
    // Use a series expansion of the BesselI[0,x] appearing in the radial integral.
    double lastSum,sum(1-std::exp(-ssSqby2)),prod(1);
    int k(0);
    while(++k < 1000) {
        lastSum = sum;
        prod *= tSqby2;
        double term(prod*boost::math::gamma_p(k+1,ssSqby2)/boost::math::tgamma(k+1));
        sum += term;
        if(std::fabs(term) < 1e-2*sum) break;
    }
    return sum*std::exp(-tSqby2);
}

int main(int argc, char **argv) {
    
    // Configure command-line option processing
    po::options_description cli("Throughput correction");
    double psfmin, psfstep;
    int plate, psfnumsteps;
    std::string inputName, prefix;
    cli.add_options()
        ("help,h", "Prints this info and exits.")
        ("verbose", "Prints additional information.")
        ("psf-min", po::value<double>(&psfmin)->default_value(1),
            "PSF fwhm minimum.")
        ("psf-num", po::value<int>(&psfnumsteps)->default_value(10),
            "Number of PSF FHWM steps.")
        ("psf-step", po::value<double>(&psfstep)->default_value(0.1,"0.1"),
            "PSF fwhm step size.")
        ("plate,p", po::value<int>(&plate)->default_value(0),
            "Plate number to use.")
        ("prefix", po::value<std::string>(&prefix)->default_value(""),
            "tpcorr file prefix.")
        ;

    // do the command line parsing now
    po::variables_map vm;
    try {
        po::store(po::parse_command_line(argc, argv, cli), vm);
        po::notify(vm);
    }
    catch(std::exception const &e) {
        std::cerr << "Unable to parse command line options: " << e.what() << std::endl;
        return -1;
    }
    if(vm.count("help")) {
        std::cout << cli << std::endl;
        return 1;
    }
    bool verbose(vm.count("verbose"));

    if (plate == 0) {
        std::cerr << "Must specify plate number!" << std::endl;
        return -2;
    }
    boost::format fmt("%f %f %f %f %f %f");
    std::string outname("tpcorr-"+boost::lexical_cast<std::string>(plate)+".dat");
    std::ofstream out(outname.c_str());

    std::string inlabelname(prefix+boost::lexical_cast<std::string>(plate)+"-label.dat");
    std::string in4000name(prefix+boost::lexical_cast<std::string>(plate)+"-4000.dat");
    std::string in5400name(prefix+boost::lexical_cast<std::string>(plate)+"-5400.dat");

    std::cout << "reading " << inlabelname << std::endl;
    std::ifstream inlabels(inlabelname.c_str());
    std::ifstream in4000(in4000name.c_str());
    std::ifstream in5400(in5400name.c_str());
    double x,y,ha,lambda;
    double dx,dy,dr4000,dr5400;
    double tp;
    int count(0);
    double psffwhm;
    while(inlabels.good() && !inlabels.eof()) {
        inlabels >> x >> y >> ha >> lambda;
        in4000 >> dx >> dy;
        dr4000 = std::sqrt(dx*dx+dy*dy)/PLATESCALE;
        in5400 >> dx >> dy;
        dr5400 = std::sqrt(dx*dx+dy*dy)/PLATESCALE;
        if(!inlabels.good() || inlabels.eof()) break;
        for(int i = 0; i < psfnumsteps; ++i) {
            psffwhm = psfmin + i*psfstep;
            if(count==0) std::cout << psffwhm << std::endl;
            tp = fiberFraction(psffwhm,dr4000)/fiberFraction(psffwhm,dr5400);
            out << fmt % x % y % ha % lambda % psffwhm % tp << std::endl;
        }
        count++;
    }
    if(verbose) {
        std::cout << "Read " << count << " lines." << std::endl;
    }
    inlabels.close();
    in4000.close();
    in5400.close();
    out.close();
    return 0;
}