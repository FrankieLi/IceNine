//
//  XDM_Monge -- a customized Monge form fitter
//
//
//

#ifndef XDM_MONGE_
#define XDM_MONGE_

#include <CGAL/Monge_via_jet_fitting.h>
#include <iostream>
CGAL_BEGIN_NAMESPACE
template< class K >
class XDM_Monge_via_jet_fitting : public Monge_via_jet_fitting< K >
{
public:

  typedef typename Monge_via_jet_fitting<K>:: Data_kernel Data_kernel;
  typedef typename Monge_via_jet_fitting<K>:: Local_kernel Local_kernel;
  
  typedef typename Monge_via_jet_fitting<K>::L2D_NTconverter L2D_NTconverter;
  typedef typename Local_kernel::FT       FT;
  typedef typename Local_kernel::Point_3  Point_3;
  typedef typename Local_kernel::Vector_3 Vector_3;
  typedef typename Monge_via_jet_fitting<K>::Aff_transformation Aff_transformation;

  typedef typename Data_kernel::FT       DFT;

  typedef typename Monge_via_jet_fitting<K>::LAVector LAVector;
  typedef typename Monge_via_jet_fitting<K>::LAMatrix LAMatrix;

  
  // overriding operator()
  typedef typename Monge_via_jet_fitting<K>::Monge_form Monge_form;

  template <class InputIterator>
  Monge_form
  operator()(InputIterator begin, InputIterator end, 
             size_t d, size_t dprime)
  {
    // precondition: on the degrees, jet and monge
    CGAL_precondition( (d >=1) && (dprime >= 1) 
                       && (dprime <= 4) && (dprime <= d) );
    this->deg = d;
    this->deg_monge = dprime;
    this->nb_d_jet_coeff = (d+1)*(d+2)/2;
    this->nb_input_pts = end - begin;
    // precondition: solvable ls system
    CGAL_precondition( nb_input_pts >= nb_d_jet_coeff );

    //Initialize
    Monge_form monge_form;
    monge_form.set_up(dprime);
    //for the system MA=Z
    LAMatrix M( this->nb_input_pts, this->nb_d_jet_coeff);
    LAVector Z( this->nb_input_pts);

    compute_PCA(begin, end);
    fill_matrix(begin, end, d, M, Z);//with precond

    if( M.number_of_rows() >=  M.number_of_columns() )
    {
      solve_linear_system(M, Z);  //correct with precond
      compute_Monge_basis(Z.vector(), monge_form);
      if ( dprime >= 3) compute_Monge_coefficients(Z.vector(), dprime, monge_form);
    }
    return monge_form;
  }
  
};

CGAL_END_NAMESPACE
#endif
