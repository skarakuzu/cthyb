
#pragma once

#include <triqs/mc_tools.hpp>
#include "./configuration.hpp"

namespace triqs_cthyb {


  struct wl_data{
    
      std::deque<std::pair<time_pt, op_desc>> opworm, opworm_dag; 
      //std::vector<std::pair<time_pt, op_desc>> opworm, opworm_dag; 
      mc_tools::random_generator &rng; 
      solve_parameters_t& params;
      int block_num;
      double cutoff = 10; // %5 tolerance
 
    enum list_space : int  { Z, G10, G11, G2 };
    struct worm_space {

        worm_space() : Nvisit(0), Vj(1.) {}
	worm_space(int Nvisit, double Vj) : Nvisit(Nvisit), Vj(Vj) {}
   	int Nvisit;
	double Vj;
    };
    
    std::map<list_space, worm_space> WS; 
    list_space current_space;


      wl_data(solve_parameters_t & params, mc_tools::random_generator &rng, int block_num) : params(params), rng(rng), block_num(block_num) {
      
      WS.insert({Z,worm_space(0,1)});
      if(params.measure_G_tau) WS.insert({G10,worm_space(0,1)});
      if(params.measure_G_tau && block_num>1) WS.insert({G11,worm_space(0,1)});
      if(params.measure_G2_tau) WS.insert({G2,worm_space(0,1)});


      //std::cout<<"Size of wormy space: "<<WS.size()<<" "<<WS.contains(G2)<<" "<<WS[G2].Nvisit<<" "<<WS[G2].Vj<<std::endl;

      }
  
    bool no_worm() { return opworm.empty() && opworm_dag.empty(); }
    
    bool chose_insert_worm() 
    {
     
    if (no_worm()) return true;
    else if(WS.contains(G2) && size_worm()==1) 
    {
     std::cout<<"Ever Here!!! "<<WS.contains(G2)<<std::endl;
     int rndn = rng(2);
     if(rndn==1) return true;
     else return false;

    }
    else return false;
	
    
//    int rndn = rng(2);
//    return rndn==1 ? 1 : 0;

    }

    void insert_worm(time_pt tau, op_desc op) { opworm.push_back(std::make_pair(tau,op)); }
    void insert_worm_dag(time_pt tau, op_desc op) { opworm_dag.push_back(std::make_pair(tau,op)); }
    //void erase_worm(int index) { opworm.erase(opworm.begin()+index); }
    //void erase_worm_dag(int index) { opworm_dag.erase(opworm_dag.begin()+index); }
    void erase_worm(int index) { opworm.pop_back(); }
    void erase_worm_dag(int index) { opworm_dag.pop_back(); }
	
    time_pt get_time(int index){ return opworm[index].first;}
    time_pt get_time_dag(int index){ return opworm_dag[index].first;}
     
    op_desc get_op(int index){ return opworm[index].second;}
    op_desc get_op_dag(int index){ return opworm_dag[index].second;}
     
    int get_block_index(int index){ return opworm[index].second.block_index;}
    int get_block_index_dag(int index){ return opworm_dag[index].second.block_index;}

    int size_worm(){return opworm.size();}
    int size_worm_dag(){return opworm_dag.size();}


    void update_current_space(){
    
    if(no_worm()) current_space = Z;
    else if(size_worm()==1 && size_worm_dag()==1) 
    {
      current_space = get_block_index(0)==0 ? G10 : G11; 
      //current_space = G10; 
    }

    else if(size_worm()==2 && size_worm_dag()==2) current_space = G2;
   
    //std::cout<<"Current space: "<<current_space<<" "<<WS[current_space].Vj<<std::endl;

    }

    void update_mu_space(){
    
    if(params.wang_landau_cycle)
    {
    WS[current_space].Vj *= params.wang_landau_lambda; 
    //WS[current_space].Nvisit += 1; 
    }
    
    double save_Vj = WS[Z].Vj;
    //for (auto it = WS.begin(); it != WS.end(); ++it) it->second.Vj /= save_Vj; 
    //for (auto &x : WS) x.second.Vj /= save_Vj;
    for (auto &[space, weight] : WS) weight.Vj /= save_Vj;
    }

    void update_visit_space(){
    
    if(params.wang_landau_cycle)
    {
    WS[current_space].Nvisit += 1; 
    }
    } 

    void reset_space_visit(){
    double save_Vj = WS[Z].Vj;
    /*
    for (auto it = WS.begin(); it != WS.end(); ++it)
    {
     it->second.Nvisit = 0; 
     it->second.Vj /= save_Vj; 
    }	
    */

    for (auto &[space, weight] : WS) { weight.Nvisit = 0; weight.Vj /= save_Vj;}
    
    //params.wang_landau_lambda = sqrt(params.wang_landau_lambda);

    }
    
    void print_space_visit(){
    for (auto it = WS.begin(); it != WS.end(); ++it)
    {
      std::cout<<"Space "<<it->first<<" has been visited "<<it->second.Nvisit<<" times with eta:"<<it->second.Vj<<std::endl;; 
    }	
    }

    double get_mu_space() {
      if(params.wang_landau_cycle) return WS[Z].Vj/WS[current_space].Vj;
      else return 1.;
    }
    
    double get_mu_G_space(int index) {
      auto G = G10; 
      //auto G = index==0 ? G10 : G11; 
      if(params.wang_landau_cycle) return WS[Z].Vj/WS[G].Vj;
      else return 1.;
    }

    bool stop_cycle()
    {
    bool res;
    
    for (auto it = WS.begin(); it != WS.end(); ++it)
    {
    res = false;
    if((abs(WS[Z].Nvisit - it->second.Nvisit) < WS[Z].Nvisit*cutoff/100. ) && params.wang_landau_lambda <=1.000000001) 
    {
     res = true; 
     params.wang_landau_lambda = 1.0;
    }
    }	
    

//    return std::all_of(WS.begin(), WS.end(), [visit_z=WS[Z].Nvisit, &cutoff, &params, this] (auto & x) {return ((abs(visit_z - x.second.Nvisit) < visit_z*cutoff/100. ) && params.wang_landau_lambda <=1); });

    return res;
    }



    };

} // namespace triqs_cthyb
