#include "../interface/kinematicTools.h"

//include c++ library function
#include <algorithm>
#include <cmath>

/*
 * compute angular separations in eta/phi units given the object directions
 */

double kinematics::deltaPhi(const double phiLeft, const double phiRight){
    static const double pi = 3.14159265358979323846;
    double phiDifference = fabs(phiLeft - phiRight);
    return std::min( phiDifference, 2*pi - phiDifference);
}

double kinematics::deltaEta(const double etaLeft, const double etaRight){
    return fabs(etaLeft - etaRight);
}

double kinematics::deltaR(const double phiLeft, const double etaLeft, const double phiRight, const double etaRight){
    double dPhi = deltaPhi(phiLeft, phiRight);
    double dEta = deltaEta(etaLeft, etaRight);
    return sqrt( dPhi*dPhi + dEta*dEta ); 
}

/*
 *find the minima and maxima for certain variables in two separate arrays of lorentz vectors
 */
double kinematics::extremum(const TLorentzVector* vecLeft, const std::vector<unsigned>& indicesLeft, const TLorentzVector* vecRight, const std::vector<unsigned>& indicesRight, 
        double (&computeVar)(const TLorentzVector&, const TLorentzVector&), const double& (& getExtremum) (const double&, const double&), double initVar)
{
    double extremum = initVar;
    //check all combinations of lorentzVectors for extrema
    for(unsigned indexLeft: indicesLeft){
        for(unsigned indexRight: indicesRight){
            double var = computeVar(vecLeft[indexLeft], vecRight[indexRight]);
            extremum = getExtremum(var, extremum);
        }
    }
    return ( (extremum != initVar) ? extremum : 0 );
}

double kinematics::minVar(const TLorentzVector* vecLeft, const std::vector<unsigned>& indicesLeft, const TLorentzVector* vecRight, const std::vector<unsigned>& indicesRight,
        double (&computeVar)(const TLorentzVector&, const TLorentzVector&) )
{
    return extremum(vecLeft, indicesLeft, vecRight, indicesRight, computeVar, std::min<double> , 99999.);
}

double kinematics::maxVar(const TLorentzVector* vecLeft, const std::vector<unsigned>& indicesLeft, const TLorentzVector* vecRight, const std::vector<unsigned>& indicesRight,
        double (&computeVar)(const TLorentzVector&, const TLorentzVector&) )
{
    return extremum(vecLeft, indicesLeft, vecRight, indicesRight, computeVar, std::max<double>, 0.);
}

/*
 * find the minima and maxima for certain variables in a single array of lorentz vectors
 */
double kinematics::extremum(const TLorentzVector* vec, const std::vector<unsigned>& indices, 
        double (&computeVar)(const TLorentzVector&, const TLorentzVector&), const double& (& getExtremum) (const double&, const double&), double initVar)
{
    double extremum = initVar;
    if(indices.size() != 0){    //safety, second loop diverges for 0 element array without this condition
        //check all combinations of different lorentzVectors for extrema
        for(unsigned outer = 0; outer < indices.size() - 1; ++outer){
            for(unsigned inner = outer + 1; inner < indices.size(); ++inner){
                double var = computeVar(vec[indices[inner]], vec[indices[outer]]);
                extremum = getExtremum(var, extremum);
            }
        }
    }
    return ( (extremum != initVar) ? extremum : 0 );
}

double kinematics::minVar(const TLorentzVector* vec, const std::vector<unsigned>& indices, double (&computeVar)(const TLorentzVector&, const TLorentzVector&) ){
    return extremum(vec, indices, computeVar, std::min<double>, 99999.);
}

double kinematics::maxVar(const TLorentzVector* vec, const std::vector<unsigned>& indices, double (&computeVar)(const TLorentzVector&, const TLorentzVector&) ){
    return extremum(vec, indices, computeVar, std::max<double>, 0.);
}


/*
 * helper functions computing several quantities given 2 TLorentzVectors
 */

double kinematics::mass(const TLorentzVector& lhs, const TLorentzVector& rhs){
    return (lhs + rhs).M();
}

double kinematics::mt(const TLorentzVector& lhs, const TLorentzVector& rhs){
    //this definition assumed the ultrarelativistic limit!
    return sqrt(2*lhs.Pt()*rhs.Pt()*( 1 - cos( lhs.Phi()-rhs.Phi() ) ) ); 
}

double kinematics::deltaR(const TLorentzVector& lhs, const TLorentzVector& rhs){
    return lhs.DeltaR(rhs);
}

double kinematics::deltaPhi(const TLorentzVector& lhs, const TLorentzVector& rhs){
    return fabs(lhs.DeltaPhi(rhs));
}

double kinematics::pt(const TLorentzVector& lhs,  const TLorentzVector& rhs){
    return (lhs + rhs).Pt();
}

/*
 * find the minima and maxima of several kinematic quantities for two separate arrays of lorentz vectors
 */
double kinematics::minMass(const TLorentzVector* vecLeft, const std::vector<unsigned>& indicesLeft, const TLorentzVector* vecRight, const std::vector<unsigned>& indicesRight){
    return minVar(vecLeft, indicesLeft, vecRight, indicesRight, mass);
}

double kinematics::maxMass(const TLorentzVector* vecLeft, const std::vector<unsigned>& indicesLeft, const TLorentzVector* vecRight, const std::vector<unsigned>& indicesRight){
    return maxVar(vecLeft, indicesLeft, vecRight, indicesRight, mass);
}

double kinematics::minDeltaR(const TLorentzVector* vecLeft, const std::vector<unsigned>& indicesLeft, const TLorentzVector* vecRight, const std::vector<unsigned>& indicesRight){
    return minVar(vecLeft, indicesLeft, vecRight, indicesRight, deltaR);
}

double kinematics::maxDeltaR(const TLorentzVector* vecLeft, const std::vector<unsigned>& indicesLeft, const TLorentzVector* vecRight, const std::vector<unsigned>& indicesRight){
    return maxVar(vecLeft, indicesLeft, vecRight, indicesRight, deltaR);
}

double kinematics::minDeltaPhi(const TLorentzVector* vecLeft, const std::vector<unsigned>& indicesLeft, const TLorentzVector* vecRight, const std::vector<unsigned>& indicesRight){
    return minVar(vecLeft, indicesLeft, vecRight, indicesRight, deltaPhi);
}

double kinematics::maxDeltaPhi(const TLorentzVector* vecLeft, const std::vector<unsigned>& indicesLeft, const TLorentzVector* vecRight, const std::vector<unsigned>& indicesRight){
    return maxVar(vecLeft, indicesLeft, vecRight, indicesRight, deltaPhi);
}

double kinematics::minMT(const TLorentzVector* vecLeft, const std::vector<unsigned>& indicesLeft, const TLorentzVector* vecRight, const std::vector<unsigned>& indicesRight){
    return minVar(vecLeft, indicesLeft, vecRight, indicesRight, mt);
}

double kinematics::maxMT(const TLorentzVector* vecLeft, const std::vector<unsigned>& indicesLeft, const TLorentzVector* vecRight, const std::vector<unsigned>& indicesRight){
    return maxVar(vecLeft, indicesLeft, vecRight, indicesRight, mt);
}

double kinematics::minPT(const TLorentzVector* vecLeft, const std::vector<unsigned>& indicesLeft, const TLorentzVector* vecRight, const std::vector<unsigned>& indicesRight){
    return minVar(vecLeft, indicesLeft, vecRight, indicesRight, pt);
}

double kinematics::maxPT(const TLorentzVector* vecLeft, const std::vector<unsigned>& indicesLeft, const TLorentzVector* vecRight, const std::vector<unsigned>& indicesRight){
    return maxVar(vecLeft, indicesLeft, vecRight, indicesRight, pt);
}

/*
 * find the minima and maxima of several kinematic quantities in one array of lorentz vectors 
 */
double kinematics::minMass(const TLorentzVector* vec, const std::vector<unsigned>& indices){
    return minVar(vec, indices, mass);
}

double kinematics::maxMass(const TLorentzVector* vec, const std::vector<unsigned>& indices){
    return maxVar(vec, indices, mass);
}

double kinematics::minDeltaR(const TLorentzVector* vec, const std::vector<unsigned>& indices){
    return minVar(vec, indices, deltaR);
}

double kinematics::maxDeltaR(const TLorentzVector* vec, const std::vector<unsigned>& indices){
    return maxVar(vec, indices, deltaR);
}

double kinematics::minDeltaPhi(const TLorentzVector* vec, const std::vector<unsigned>& indices){
    return minVar(vec, indices, deltaPhi);
}

double kinematics::maxDeltaPhi(const TLorentzVector* vec, const std::vector<unsigned>& indices){
    return maxVar(vec, indices, deltaPhi);
}

double kinematics::minMT(const TLorentzVector* vec, const std::vector<unsigned>& indices){
    return minVar(vec, indices, mt);
}

double kinematics::maxMT(const TLorentzVector* vec, const std::vector<unsigned>& indices){
    return maxVar(vec, indices, mt);
}

double kinematics::minPT(const TLorentzVector* vec, const std::vector<unsigned>& indices){
    return minVar(vec, indices, pt);
}

double kinematics::maxPT(const TLorentzVector* vec, const std::vector<unsigned>& indices){
    return maxVar(vec, indices, pt);
}
