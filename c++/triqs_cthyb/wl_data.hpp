
#pragma once

#include "./configuration.hpp"

namespace triqs_cthyb {


  struct wl_data{
    
      //using opworm_t = std::map<time_pt, op_desc>;
      std::vector<std::pair<time_pt, op_desc>> opworm; 
      std::vector<std::pair<time_pt, op_desc>> opworm_dag; 
      double wl_lambda;
      std::vector<double> V_l;
      //opworm_t opworm;
    
    
      wl_data(double lmbda) : wl_lambda(lmbda) {}
   /* 
    void insert_worm(time_pt tau, op_desc op) { opworm.insert({tau, op}); }
    void replace_worm(time_pt tau, op_desc op) { opworm[tau] = op; }
    void erase_worm(time_pt const &t) { opworm.erase(t); }
    void clear_worm() { opworm.clear(); }
    int number_of_worm() { return opworm.size()/2; }
    */
    bool no_worm() { return opworm.empty() && opworm_dag.empty(); }
    void insert_worm(time_pt tau, op_desc op) { opworm.push_back(std::make_pair(tau,op)); }
    void insert_worm_dag(time_pt tau, op_desc op) { opworm_dag.push_back(std::make_pair(tau,op)); }
    void erase_worm(int index) { opworm.erase(opworm.begin()+index); }
    void erase_worm_dag(int index) { opworm_dag.erase(opworm_dag.begin()+index); }
	
    time_pt get_time(int index){ return opworm[0].first;}
    time_pt get_time_dag(int index){ return opworm_dag[0].first;}
     
    int get_block_index(int index){ return opworm[0].second.block_index;}
    int get_block_index_dag(int index){ return opworm_dag[0].second.block_index;}

    int size_worm(){return opworm.size();}
    int size_worm_dag(){return opworm_dag.size();}

    };

} // namespace triqs_cthyb
