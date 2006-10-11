#ifndef DF_REPORT_COLLECTION_HPP
#define DF_REPORT_COLLECTION_HPP
#include "port.h"

#include <vector>
using namespace std;

#include "Util_Abstract_Group.hpp"
#include "DF_Abstract_Report.hpp"

namespace DFLib
{
    class ReportCollection : public DFLib::Abstract::Group
    {
    private: 
        vector<DFLib::Abstract::Report * > theReports;
        vector<double> evaluationPoint;
        bool f_is_valid;
        bool g_is_valid;
        bool h_is_valid;
        double function_value;
        vector<double> gradient;
        vector<vector<double> > hessian;

        // Declare the copy constructor and assignment operators, but
        // don't define them.  We should *never* copy a collection or attempt
        // to assign one to another.  This makes it illegal to do so.
        ReportCollection(ReportCollection &right);
        ReportCollection &operator=(ReportCollection &right);

    public:
        CPL_DLL ReportCollection();

        /// \brief DF Report Collection destructor
        ///
        ///  Never destroys the objects in its vector, as they might be getting
        ///  used for something else by caller
        CPL_DLL ~ReportCollection();

        /// \brief destroy all reports stored in collection
        ///
        /// Provided in case our caller does NOT need the stored pointers for
        /// something else, and doesn't want to keep track of them.
        CPL_DLL void deleteReports();

        /// \brief Add a DF report to the collection
        ///
        /// \return this report's number in the collection.
        CPL_DLL int addReport(DFLib::Abstract::Report * aReport);

        /// \brief return the fix cut average of this collection's reports
        ///
        /// A fix cut is the intersection of two DF reports.  The Fix Cut 
        /// Average is the average of all fix cuts from all pairs of reports 
        /// in the collection.  Fix cuts at shallow angles can be excluded by 
        /// specifying a non-zero value for minAngle (in degrees).
        ///
        /// \param FCA Returned fix cut average
        /// \param FCA_stddev standard deviation of fix cuts
        /// \param minAngle reports whose fix cut occur at less than this angle will not be included in the average.
        CPL_DLL bool computeFixCutAverage(vector<double> &FCA, 
                                  vector<double> &FCA_stddev,
                                  double minAngle=0);

        /*!
          \brief Computes least squares solution of DF problem.
       
         The least squares solution is the point that has the minimum
         orthogonal distance to all bearing lines.
        
         The least squares solution is given by

         \f$ 
         P_{LS} = (A^TA)^{-1}A^Tb 
         \f$

         where  \f$A\f$ is the \f$Nx2\f$ matrix whose \f$k^{th}\f$ row
         is the unit vector orthogonal to the bearing line from receiver
         \f$k\f$.  Row \f$k\f$ of A is therefore the vector
         \f$[\cos(\theta_k), -\sin(\theta_k)]\f$ where \f$\theta_k\f$
         is the bearing from station \f$k\f$ measured clockwise from North.
         The \f$T\f$ superscript denotes matrix transpose.
        
         \f$b\f$ is the \f$N\f$ element vector whose \f$k^{th}\f$ element
         is the projection of the position vector of the \f$k^{th}\f$ 
         receiver onto the unit vector orthogonal to that receiver's 
         bearing line:  
         \f$
           b_k = [\cos(\theta_k),-sin(\theta_k)]\cdot[r_{kx},r_{ky}]
         \f$
         where \f${\bf r}_k\f$ is the position vector of the \f$k^{th}\f$
         receiver.
        */
        CPL_DLL void computeLeastSquaresFix(vector<double> &LS_Fix);
        CPL_DLL double computeCostFunction(vector<double> &evaluationPoint);
        CPL_DLL void computeCostFunctionAndGradient(vector<double> &evaluationPoint,
                                            double &f,
                                            vector<double> &gradf);
        CPL_DLL void computeCostFunctionAndHessian(vector<double> &evaluationPoint,
                                                 double &f,
                                                 vector<double> &gradf,
                                                 vector<vector<double> > &h);

        inline CPL_DLL int size() {return theReports.size();};
        inline CPL_DLL const DFLib::Abstract::Report * const getReport(int i) 
        { return (theReports[i]); };

        inline CPL_DLL const vector<double> &getReceiverLocationXY(int i) 
        { return (theReports[i]->getReceiverLocation()); };

        inline virtual CPL_DLL void setEvaluationPoint(vector<double> &ep)
        {
            evaluationPoint = ep;
            f_is_valid=false;
            g_is_valid=false;
            h_is_valid=false;
        }

        /// \return function value
        inline virtual CPL_DLL double getFunctionValue()
        {
            if (!f_is_valid)
            {
                function_value=computeCostFunction(evaluationPoint);
                f_is_valid=true;
            }
            return function_value;
        }

        /// \param g gradient returned
        /// \return function value
        inline virtual CPL_DLL double getFunctionValueAndGradient(vector<double> &g)
        {
            if (!g_is_valid)
            {
                computeCostFunctionAndGradient(evaluationPoint,
                                               function_value,gradient);
                f_is_valid=true;
                g_is_valid=true;
            }
            g=gradient;
            return function_value;
        }

        /// \param g gradient returned
        /// \param h hessian returned
        /// \return function value
        inline virtual CPL_DLL double 
        getFunctionValueAndHessian(vector<double> &g,
                                   vector<vector<double> > &h)
        {
            if (!h_is_valid)
                {
                    computeCostFunctionAndHessian(evaluationPoint,
                                                  function_value,gradient,
                                                  hessian);
                    f_is_valid=true;
                    g_is_valid=true;
                    h_is_valid=true;
                }
            g=gradient;
            h=hessian;
            return function_value;
        }
    };
}
#endif // DF_REPORT_COLLECTION_HPP
