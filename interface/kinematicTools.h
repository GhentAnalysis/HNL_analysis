/*
 * functions to do some common kinematic computations
 */


//include ROOT classes
#include "TLorentzVector.h"


namespace kinematics{
    /*
     * functions to compute angular separations in eta/phi units given the object directions
     */
    double deltaPhi(const double, const double);
    double deltaEta(const double, const double); 
    double deltaR(const double, const double, const double, const double);

    /*
     * functions computing several quantities given 2 TLorentzVectors
     */

    double mass(const TLorentzVector&, const TLorentzVector&);
    double mt(const TLorentzVector&, const TLorentzVector&);
    double deltaR(const TLorentzVector&, const TLorentzVector&);
    double deltaPhi(const TLorentzVector&, const TLorentzVector&);
    double pt(const TLorentzVector&, const TLorentzVector&);

    /*
     * functions to compute kinematic extrema in the event
     */

    double minMass_OS(const std::vector<TLorentzVector>&, const std::vector<int>&, const std::vector<unsigned>&, int vector_pair[2]);
    double minVar_OS(const std::vector<TLorentzVector>&, const std::vector<int>&, const std::vector<unsigned>&, int vector_pair[2],double (&computeVar)(const TLorentzVector&, const TLorentzVector&));
    double extremum_OS(const std::vector<TLorentzVector>& , const std::vector<int>& , const std::vector<unsigned>& , int vector_pair[2], double (&computeVar)(const TLorentzVector&, const TLorentzVector&), const double& (& getExtremum) (const double&, const double&), double initVar);

    
    //parameters are 2 separate arrays of TLorentzVectors, and a vector with the relevant indices to access the correct elements
    double minMass(const TLorentzVector*, const std::vector<unsigned>&, const TLorentzVector*, const std::vector<unsigned>&);
    double maxMass(const TLorentzVector*, const std::vector<unsigned>&, const TLorentzVector*, const std::vector<unsigned>&);
    double minDeltaR(const TLorentzVector*, const std::vector<unsigned>&, const TLorentzVector*, const std::vector<unsigned>&);
    double maxDeltaR(const TLorentzVector*, const std::vector<unsigned>&, const TLorentzVector*, const std::vector<unsigned>&);
    double minDeltaPhi(const TLorentzVector*, const std::vector<unsigned>&, const TLorentzVector*, const std::vector<unsigned>&);
    double maxDeltaPhi(const TLorentzVector*, const std::vector<unsigned>&, const TLorentzVector*, const std::vector<unsigned>&);
    double minMT(const TLorentzVector*, const std::vector<unsigned>&, const TLorentzVector*, const std::vector<unsigned>&);
    double maxMT(const TLorentzVector*, const std::vector<unsigned>&, const TLorentzVector*, const std::vector<unsigned>&);
    double minPT(const TLorentzVector*, const std::vector<unsigned>&, const TLorentzVector*, const std::vector<unsigned>&);
    double maxPT(const TLorentzVector*, const std::vector<unsigned>&, const TLorentzVector*, const std::vector<unsigned>&);

    //parameters are one TLorentzVector, and a vector with the relevant indices to access the correct elements
    double minMass(const TLorentzVector*, const std::vector<unsigned>&);
    double maxMass(const TLorentzVector*, const std::vector<unsigned>&);
    double minDeltaR(const TLorentzVector*, const std::vector<unsigned>&);
    double maxDeltaR(const TLorentzVector*, const std::vector<unsigned>&);
    double minDeltaPhi(const TLorentzVector*, const std::vector<unsigned>&);
    double maxDeltaPhi(const TLorentzVector*, const std::vector<unsigned>&);
    double minMT(const TLorentzVector*, const std::vector<unsigned>&);
    double maxMT(const TLorentzVector*, const std::vector<unsigned>&);
    double minPT(const TLorentzVector*, const std::vector<unsigned>&);
    double maxPT(const TLorentzVector*, const std::vector<unsigned>&);

    //helper functions to avoid duplicate loop code
    /*
     * find the minima and maxima for certain variables in two separate arrays of lorentz vectors
     */
    double extremum(const TLorentzVector*, const std::vector<unsigned>&, const TLorentzVector*, const std::vector<unsigned>&,
                    double (&computeVar)(const TLorentzVector&, const TLorentzVector&), const double& (& getExtremum) (const double&, const double&), double initVar);
    double minVar(const TLorentzVector*, const std::vector<unsigned>&, const TLorentzVector*, const std::vector<unsigned>&,
                  double (&computeVar)(const TLorentzVector&, const TLorentzVector&) );
    double maxVar(const TLorentzVector*, const std::vector<unsigned>&, const TLorentzVector*, const std::vector<unsigned>&,
                  double (&computeVar)(const TLorentzVector&, const TLorentzVector&) );

    /*
     * find the minima and maxima for certain variables in a single array of lorentz vectors
     */ 
    double extremum(const TLorentzVector*, const std::vector<unsigned>&,
                    double (&computeVar)(const TLorentzVector&, const TLorentzVector&), const double& (& getExtremum) (const double&, const double&), double initVar);
    double minVar(const TLorentzVector*, const std::vector<unsigned>&, double (&computeVar)(const TLorentzVector&, const TLorentzVector&) );
    double maxVar(const TLorentzVector*, const std::vector<unsigned>&, double (&computeVar)(const TLorentzVector&, const TLorentzVector&) );

}   
