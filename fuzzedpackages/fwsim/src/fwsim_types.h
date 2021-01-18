#ifndef _haptools_FWSIM_H
#define _haptools_FWSIM_H

#include <Rcpp.h>

// [[Rcpp::plugins(cpp11)]]

class Individual {
  private:
    int m_id;
    int m_individual;
    int m_generation;
    int m_founder_id;
    std::vector<int> m_haplotype;
    Individual* m_parent;
    std::vector<Individual*> m_children;

  public:
    static void cleanup_lineage(Individual* node) {
      Individual* parent = node->get_parent();
      
      if (node->get_children().size() == 0) {
        Individual* parent = node->get_parent();
        
        if (parent) {
          parent->remove_child(node);     
          
          // Put here to avoid deleting founders
          //Rcout << "deletes " << node->get_label() << std::endl;
          delete node;
          
          Individual::cleanup_lineage(parent);
        }
      }
    }

    Individual(int id, int generation, int individual, Individual* parent, std::vector<int> haplotype) {
      m_id = id;
      m_generation = generation;
      m_individual = individual;
      
      if (parent != NULL) {
        m_founder_id = parent->get_founder_id();
      } else {
        m_founder_id = individual;
      }
      
      m_haplotype = haplotype;
      m_parent = parent;
    }
    
    ~Individual () {
      //m_parent->remove_child(this);
    }
    
    /*
    bool operator ==(Individual& i2) const {
      Rcout << "ASD" << std::endl;
      return i2.get_id() == m_id;
    }
    */
    
    void cleanup_haplotype() {
      m_haplotype.clear();
    }    
    
    int get_founder_id() {
      return m_founder_id;
    }
    
    void add_child(Individual* child) {
      m_children.push_back(child);
    }
    
    void remove_child(Individual* child) {
      if (m_children.size() == 0) {
        return;
      } 

      int id = child->get_id();
      
      for (std::vector<Individual*>::iterator iter = m_children.begin(); iter != m_children.end(); ++iter ) {
        if ((*iter)->get_id() == id) {
          //Rcout << "Removes id = " << id << std::endl;
          m_children.erase(iter);
          break;
        }
      }
    }
    
    std::vector<int> get_haplotype() const {
      return m_haplotype;
    }
    
    Individual* get_parent() const {
      return m_parent;
    }
    
    std::vector<Individual*> get_children() const {
      return m_children;      
    }
    
    int get_id() const {
      return m_id;
    }
    
    int get_generation() const {
      return m_generation;
    }
    
    std::string get_label() const { 
      std::ostringstream oss;
      
      if (m_haplotype.size() > 0) {
        oss << ": (";
        std::copy(m_haplotype.begin(), m_haplotype.end() - 1, std::ostream_iterator<int>(oss, ","));
        oss << m_haplotype.back();
        oss << ")";
      }
      
      std::ostringstream ss;
      //ss << m_generation << "_" << m_individual;
      ss << "[" << m_id << "/" << m_founder_id << "] " << m_generation << "_" << m_individual << "" << oss.str();
      std::string s = ss.str();
      return s;
    }
};


class SimulatedGenealogy {
  private:
    int m_pop_size;
    std::vector<Individual*> m_population;
    std::vector<Individual*> m_init_population;

  public:
    SimulatedGenealogy(std::vector<Individual*> population, std::vector<Individual*> init_population) {
      m_pop_size = population.size();

      if (m_pop_size == 0) {
        throw std::invalid_argument("population must have size of at least 1");
      }
      
      if (init_population.size() != m_pop_size) {
        throw std::invalid_argument("population and init_population must have same size");
      }
      
      m_population = population;
      m_init_population = init_population;
    }
    
    ~SimulatedGenealogy() {
      //std::cout << "SimulatedGenealogy died" << std::endl;
      /*
      // Clean up
      for (int individual = 0; individual < m_pop_size; individual++) {
        Individual::cleanup_lineage(m_population[individual]);
      }
        
      for (int individual = 0; individual < m_pop_size; ++individual) {
        delete m_init_population[individual];    
      } 
      */      
    }
    
    int get_pop_size() const {
      return m_pop_size;
    }
    
    // give access to private data
    friend std::ostream& operator<<(std::ostream &strm, const SimulatedGenealogy &simres);
};

std::ostream& operator<<(std::ostream &strm, const SimulatedGenealogy &simres);

#endif

