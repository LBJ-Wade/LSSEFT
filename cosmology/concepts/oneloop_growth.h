//
// Created by David Seery on 17/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_ONELOOP_H
#define LSSEFT_ONELOOP_H


#include <memory>
#include <vector>

#include "database/tokens.h"
#include "database/z_database.h"

#include "boost/timer/timer.hpp"


struct oneloop_growth_record
  {
    double g;
    double A;
    double B;
    double D;
    double E;
    double F;
    double G;
    double J;
    
    double f;
    double fA;
    double fB;
    double fD;
    double fE;
    double fF;
    double fG;
    double fJ;
  };


typedef std::pair< const z_token&, oneloop_growth_record > oneloop_value;


namespace oneloop_growth_impl
  {

    template <typename RecordIterator, typename ConstRecordIterator, typename ValueIterator, typename ConstValueIterator, bool is_const_iterator=true>
    class generic_tokenized_database_iterator: public std::iterator< std::bidirectional_iterator_tag, oneloop_value >
      {

      private:

        //! type alias for raw record iterator
        typedef typename std::conditional< is_const_iterator, ConstRecordIterator, RecordIterator >::type record_iterator;

        //! type alias for raw value iterator
        typedef typename std::conditional< is_const_iterator, ConstValueIterator, ValueIterator >::type value_iterator;


        // CONSTRUCTOR, DESTRUCTOR

      public:

        //! default constructor: points to nothing when it is constructed
        generic_tokenized_database_iterator() = default;

        //! value constructor; points to a given element in the transfer function sample
        generic_tokenized_database_iterator(record_iterator r,
                                            value_iterator g, value_iterator A, value_iterator B, value_iterator D,
                                            value_iterator E, value_iterator F, value_iterator G, value_iterator J,
                                            value_iterator f, value_iterator fA, value_iterator fB, value_iterator fD,
                                            value_iterator fE, value_iterator fF, value_iterator fG, value_iterator fJ)
          : record_iter(r),
            g_iter(g),
            A_iter(A),
            B_iter(B),
            D_iter(D),
            E_iter(E),
            F_iter(F),
            G_iter(G),
            J_iter(J),
            f_iter(f),
            fA_iter(fA),
            fB_iter(fB),
            fD_iter(fD),
            fE_iter(fE),
            fF_iter(fF),
            fG_iter(fG),
            fJ_iter(fJ)
          {
          }

        //! copy constructor; allows implicit conversion from a regular iterator to a const iterator
        generic_tokenized_database_iterator(const generic_tokenized_database_iterator<RecordIterator, ConstRecordIterator, ValueIterator, ConstValueIterator, false>& obj)
          : record_iter(obj.record_iter),
            g_iter(obj.g_iter),
            A_iter(obj.A_iter),
            B_iter(obj.B_iter),
            D_iter(obj.D_iter),
            E_iter(obj.E_iter),
            F_iter(obj.F_iter),
            G_iter(obj.G_iter),
            J_iter(obj.J_iter),
            f_iter(obj.f_iter),
            fA_iter(obj.fA_iter),
            fB_iter(obj.fB_iter),
            fD_iter(obj.fD_iter),
            fE_iter(obj.fE_iter),
            fF_iter(obj.fF_iter),
            fG_iter(obj.fG_iter),
            fJ_iter(obj.fJ_iter)
          {
          }


        // OVERLOAD COMPARISON OPERATORS

      public:

        //! equality comparison
        bool operator==(const generic_tokenized_database_iterator& obj) const
          {
            // all should be in step, so need only compare one of them
            return(this->record_iter == obj.record_iter);
          }

        //! inequality comparison
        bool operator!=(const generic_tokenized_database_iterator& obj) const
          {
            // all iterators should be in step, so need only compare one of them
            return(this->record_iter != obj.record_iter);
          }


        // DEREFERENCING

      public:

        //! dereference iterator to get value
        oneloop_value operator*() const
          {
            oneloop_growth_record rec;

            rec.g = *this->g_iter;
            rec.A = *this->A_iter;
            rec.B = *this->B_iter;
            rec.D = *this->D_iter;
            rec.E = *this->E_iter;
            rec.F = *this->F_iter;
            rec.G = *this->G_iter;
            rec.J = *this->J_iter;
            
            rec.f = *this->f_iter;
            rec.fA = *this->fA_iter;
            rec.fB = *this->fB_iter;
            rec.fD = *this->fD_iter;
            rec.fE = *this->fE_iter;
            rec.fG = *this->fG_iter;
            rec.fF = *this->fF_iter;
            rec.fJ = *this->fJ_iter;

            return oneloop_value(this->record_iter->get_token(), rec);
          }


        // INCREMENT, DECREMENT

      public:

        //! prefix decrement
        generic_tokenized_database_iterator& operator--()
          {
            --this->g_iter;
            --this->A_iter;
            --this->B_iter;
            --this->D_iter;
            --this->E_iter;
            --this->F_iter;
            --this->G_iter;
            --this->J_iter;
            --this->record_iter;
            return(*this);
          }

        //! postfix decrement
        generic_tokenized_database_iterator& operator--(int)
          {
            const generic_tokenized_database_iterator old(*this);
            --(*this);
            return(old);
          }

        //! prefix increment
        generic_tokenized_database_iterator& operator++()
          {
            ++this->g_iter;
            ++this->A_iter;
            ++this->B_iter;
            ++this->D_iter;
            ++this->E_iter;
            ++this->F_iter;
            ++this->G_iter;
            ++this->J_iter;
            
            ++this->f_iter;
            ++this->fA_iter;
            ++this->fB_iter;
            ++this->fD_iter;
            ++this->fE_iter;
            ++this->fF_iter;
            ++this->fG_iter;
            ++this->fJ_iter;

            ++this->record_iter;

            return(*this);
          }

        //! postfix increment
        generic_tokenized_database_iterator& operator++(int)
          {
            const generic_tokenized_database_iterator old(*this);
            ++(*this);
            return(old);
          }

        // make the const version a friend of the non-const version,
        // so the copy constructor can access its private members during implicit conversion
        friend class generic_tokenized_database_iterator<RecordIterator, ConstRecordIterator, ValueIterator, ConstValueIterator, true>;


        // INTERNAL DATA

      private:

        //! iterator into redshift database
        record_iterator record_iter;

        //! iterator into linear growth factor sample
        value_iterator g_iter;

        //! iterator into A sample
        value_iterator A_iter;

        //! iterator into B sample
        value_iterator B_iter;

        //! iterator into D sample
        value_iterator D_iter;

        //! iterator into E sample
        value_iterator E_iter;

        //! iterator into F sample
        value_iterator F_iter;

        //! iterator into G sample
        value_iterator G_iter;

        //! iterator into J sample
        value_iterator J_iter;
        
        //! iterator into linear growth rate sample
        value_iterator f_iter;
        
        //! iterator into fA sample
        value_iterator fA_iter;
    
        //! iterator into fB sample
        value_iterator fB_iter;
    
        //! iterator into fD sample
        value_iterator fD_iter;
    
        //! iterator into fE sample
        value_iterator fE_iter;
    
        //! iterator into fF sample
        value_iterator fF_iter;
    
        //! iterator into fG sample
        value_iterator fG_iter;
    
        //! iterator into fJ sample
        value_iterator fJ_iter;

      };


  }   // namespace oneloop_growth_impl


class oneloop_growth
  {

    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! constructor
    oneloop_growth(z_database& z);

    //! destructor is default
    ~oneloop_growth() = default;


    // ITERATORS

  public:

    //! type alias for non-const iterator
    typedef oneloop_growth_impl::generic_tokenized_database_iterator<z_database::reverse_record_iterator, z_database::const_reverse_record_iterator,
                                                                     std::vector<double>::iterator, std::vector<double>::const_iterator, false> iterator;

    //! type alias for const iterator
    typedef oneloop_growth_impl::generic_tokenized_database_iterator<z_database::reverse_record_iterator, z_database::const_reverse_record_iterator,
                                                                     std::vector<double>::iterator, std::vector<double>::const_iterator, true> const_iterator;

    iterator begin()
      {
        return(iterator(this->z_db.record_rbegin(),
                        this->g_linear->begin(), this->A->begin(), this->B->begin(), this->D->begin(), this->E->begin(), this->F->begin(), this->G->begin(), this->J->begin(),
                        this->f_linear->begin(), this->fA->begin(), this->fB->begin(), this->fD->begin(), this->fE->begin(), this->fF->begin(), this->fG->begin(), this->fJ->begin()));
      }

    iterator end()
      {
        return(iterator(this->z_db.record_rend(),
                        this->g_linear->end(), this->A->end(), this->B->end(), this->D->end(), this->E->end(), this->F->end(), this->G->end(), this->J->end(),
                        this->f_linear->end(), this->fA->end(), this->fB->end(), this->fD->end(), this->fE->end(), this->fF->end(), this->fG->end(), this->fJ->end()));
      }

    const_iterator begin() const
      {
        return(const_iterator(this->z_db.record_crbegin(),
                              this->g_linear->cbegin(), this->A->cbegin(), this->B->cbegin(), this->D->cbegin(), this->E->cbegin(), this->F->cbegin(), this->G->cbegin(), this->J->cbegin(),
                              this->f_linear->cbegin(), this->fA->cbegin(), this->fB->cbegin(), this->fD->cbegin(), this->fE->cbegin(), this->fF->cbegin(), this->fG->cbegin(), this->fJ->cbegin()));
      }

    const_iterator end() const
      {
        return(const_iterator(this->z_db.record_crend(),
                              this->g_linear->cend(), this->A->cend(), this->B->cend(), this->D->cend(), this->E->cend(), this->F->cend(), this->G->cend(), this->J->cend(),
                              this->f_linear->cend(), this->fA->cend(), this->fB->cend(), this->fD->cend(), this->fE->cend(), this->fF->cend(), this->fG->cend(), this->fJ->cend()));
      }

    const_iterator cbegin() const
      {
        return(const_iterator(this->z_db.record_crbegin(),
                              this->g_linear->cbegin(), this->A->cbegin(), this->B->cbegin(), this->D->cbegin(), this->E->cbegin(), this->F->cbegin(), this->G->cbegin(), this->J->cbegin(),
                              this->f_linear->cbegin(), this->fA->cbegin(), this->fB->cbegin(), this->fD->cbegin(), this->fE->cbegin(), this->fF->cbegin(), this->fG->cbegin(), this->fJ->cbegin()));
      }

    const_iterator cend() const
      {
        return(const_iterator(this->z_db.record_crend(),
                              this->g_linear->cend(), this->A->cend(), this->B->cend(), this->D->cend(), this->E->cend(), this->F->cend(), this->G->cend(), this->J->cend(),
                              this->f_linear->cend(), this->fA->cend(), this->fB->cend(), this->fD->cend(), this->fE->cend(), this->fF->cend(), this->fG->cend(), this->fJ->cend()));
      }


    // INTERFACE

  public:

    //! store components
    void push_back(double g, double A, double B, double D, double E, double F, double G, double J,
                   double f, double fA, double fB, double fD, double fE, double fF, double fG, double fJ);


    // METADATA

  public:

    //! store integration time
    void set_integration_metadata(boost::timer::nanosecond_type t, size_t s);

    //! get integration time
    boost::timer::nanosecond_type get_integration_time() const { return(this->integration_time); }

    //! get number of steps used by integrator
    size_t get_integration_steps() const { return(this->steps); }


    // INTERNAL DATA

  private:

    // CONFIGURATION DATA

    //! reference to redshift database
    z_database& z_db;


    // ONE-LOOP FUNCTIONS

    // these are managed using std::unique_ptr<>s to control their lifetime

    //! linear growth factor g(z)
    std::unique_ptr< std::vector<double> > g_linear;

    //! A growth factor
    std::unique_ptr< std::vector<double> > A;

    //! B growth factor
    std::unique_ptr< std::vector<double> > B;

    //! D growth factor
    std::unique_ptr< std::vector<double> > D;

    //! E growth factor
    std::unique_ptr< std::vector<double> > E;

    //! F growth factor
    std::unique_ptr< std::vector<double> > F;

    //! G growth factor
    std::unique_ptr< std::vector<double> > G;

    //! J growth factor
    std::unique_ptr< std::vector<double> > J;

    
    //! linear growth rate f(z)
    std::unique_ptr< std::vector<double> > f_linear;
    
    //! A growth rate
    std::unique_ptr< std::vector<double> > fA;
    
    //! B growth rate
    std::unique_ptr< std::vector<double> > fB;
    
    //! D growth rate
    std::unique_ptr< std::vector<double> > fD;
    
    //! E growth rate
    std::unique_ptr< std::vector<double> > fE;
    
    //! F growth rate
    std::unique_ptr< std::vector<double> > fF;
    
    //! G growth rate
    std::unique_ptr< std::vector<double> > fG;
    
    //! J growth rate
    std::unique_ptr< std::vector<double> > fJ;

    
    // METADATA

    //! time taken to perform integration
    boost::timer::nanosecond_type integration_time;

    //! number of steps used by integrator
    size_t steps;

  };


#endif //LSSEFT_ONELOOP_H
