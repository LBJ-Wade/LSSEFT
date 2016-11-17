//
// Created by David Seery on 22/11/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_DATA_MANAGER_IMPL_H
#define LSSEFT_DATA_MANAGER_IMPL_H


#include "k_database.h"
#include "IR_cutoff_database.h"
#include "UV_cutoff_database.h"
#include "IR_resum_database.h"


namespace data_manager_impl
  {

    class loop_momentum_configuration
      {

      public:

        loop_momentum_configuration(k_database::const_record_iterator& _k,
                                    UV_cutoff_database::const_record_iterator& _UV_cutoff,
                                    IR_cutoff_database::const_record_iterator& _IR_cutoff)
          : k(_k),
            UV_cutoff(_UV_cutoff),
            IR_cutoff(_IR_cutoff)
          {
          }

        ~loop_momentum_configuration() = default;

        k_database::const_record_iterator k;
        UV_cutoff_database::const_record_iterator UV_cutoff;
        IR_cutoff_database::const_record_iterator IR_cutoff;
      };
    
    
    // specialize operator== to test for equality of loop_momentum_configuration items
    inline bool operator==(const loop_momentum_configuration& A, const loop_momentum_configuration& B)
      {
        return A.k == B.k && A.UV_cutoff == B.UV_cutoff && A.IR_cutoff == B.IR_cutoff;
      }
    
    
    class resummed_Pk_configuration
      {
      
      public:
    
        resummed_Pk_configuration(k_database::const_record_iterator& _k,
                                  UV_cutoff_database::const_record_iterator& _UV_cutoff,
                                  IR_cutoff_database::const_record_iterator& _IR_cutoff,
                                  IR_resum_database::const_record_iterator& _IR_resum)
          : k(_k),
            UV_cutoff(_UV_cutoff),
            IR_cutoff(_IR_cutoff),
            IR_resum(_IR_resum)
          {
          }
        
        ~resummed_Pk_configuration() = default;
        
        k_database::const_record_iterator k;
        UV_cutoff_database::const_record_iterator UV_cutoff;
        IR_cutoff_database::const_record_iterator IR_cutoff;
        IR_resum_database::const_record_iterator IR_resum;
      };
    
    
    // specialize operator== to test for equality of loop_momentum_configuration items
    inline bool operator==(const resummed_Pk_configuration& A, const resummed_Pk_configuration& B)
      {
        return A.k == B.k && A.UV_cutoff == B.UV_cutoff && A.IR_cutoff == B.IR_cutoff && A.IR_resum == B.IR_resum;
      }

  }


namespace std
  {
    
    namespace hash_impl
      {
    
        // base case for variadic recursion
        inline void hash_combine(std::size_t& seed)
          {
          }
    
    
        // inductive case for variadic recursion
        template <typename T, typename... Rest>
        inline void hash_combine(std::size_t& seed, const T& v, Rest... rest)
          {
            std::hash<T> hasher;
            seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            hash_combine(seed, rest...);
          }
        
      }   // namespace hash_impl


    // specialize std::hash() to loop_momentum_configuration so that it can be used with std::unordered_set
    template <>
    class hash<::data_manager_impl::loop_momentum_configuration>
      {
        
        // CONSTRUCTOR, DESTRUCTOR
        
      public:
        
        //! constructor is default
        hash() = default;
        
        //! destructor is default
        ~hash() = default;
        
        
        // IMPLEMENTATION
        
      public:
        
        //! hash operation
        size_t operator()(const ::data_manager_impl::loop_momentum_configuration& config) const
          {
            std::size_t hash = 0;
            hash_impl::hash_combine(hash, config.k->get_token(), config.UV_cutoff->get_token(), config.IR_cutoff->get_token());
            
            return hash;
          }
        
      };
    
    
    // specialize std::hash() to resummed_Pk_configuration so that it can be used with std::unordered_set
    template <>
    class hash<::data_manager_impl::resummed_Pk_configuration>
      {
        
        // CONSTRUCTOR, DESTRUCTOR
      
      public:
        
        //! constructor is default
        hash() = default;
        
        //! destructor is default
        ~hash() = default;
        
        
        // IMPLEMENTATION
      
      public:
        
        //! hash operation
        size_t operator()(const ::data_manager_impl::resummed_Pk_configuration& config) const
          {
            std::size_t hash = 0;
            hash_impl::hash_combine(hash, config.k->get_token(), config.UV_cutoff->get_token(),
                                    config.IR_cutoff->get_token(), config.IR_resum->get_token());
            
            return hash;
          }
        
      };
    
  }   // namespace std


#endif //LSSEFT_DATA_MANAGER_IMPL_H
