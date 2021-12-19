
namespace triqs_cthyb {


  struct wl_data{
    
      double wl_lambda;
      std::vector<double> V_l;
    
      wl_data(double lmbda) : wl_lambda(lmbda) {}
    };

} // namespace triqs_cthyb
