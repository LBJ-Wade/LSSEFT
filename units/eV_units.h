//
// Created by David Seery on 11/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_UNITS_H
#define LSSEFT_UNITS_H


#include "boost/math/constants/constants.hpp"
#include "boost/serialization/serialization.hpp"


namespace eV_units
  {

    // set up a template class representing a unit in a natural system where c=hbar=1
    // Here, we use eV as the base unit
    // That only leaves the reduced Planck mass Mp as a dimensionful object,
    // and we express everything else in terms of it

    // adapted from Bjarne Stroustrup's Going Native 2012 example
    // https://channel9.msdn.com/Events/GoingNative/GoingNative-2012/Keynote-Bjarne-Stroustrup-Cpp11-Style

    template <int m>
    struct unit
      {
        enum { MassDimension=m };
      };


    template <typename Unit>
    class value
      {

        // CONSTRUCTOR, DESTRUCTOR

      public:

        //! constexpr constructor means that this class can be used to build compile-time constant expressions
        constexpr value(double d)
          : val(d)
          {
          }

        ~value() = default;


        // INTERFACE

      public:

        //! allow implicit conversion to double
        explicit constexpr operator double() const { return(this->val); }


        // INTERNAL DATA
        // not made private to avoid overloads of arithmetic operators etc. becoming too heavyweight

      public:

        //! numerical value
        double val;


        // enable boost::serialization support, and hence automated packing for transmission over MPI
        friend class boost::serialization::access;

        template <typename Archive>
        void serialize(Archive& ar, unsigned int version)
          {
            ar & val;
          }

      };

    using eV_unit          = unit<1>;
    using eV2_unit         = unit<2>;
    using eV3_unit         = unit<3>;
    using inverse_eV_unit  = unit<-1>;
    using inverse_eV2_unit = unit<-2>;
    using inverse_eV3_unit = unit<-3>;

    using energy           = value<eV_unit>;
    using energy2          = value<eV2_unit>;
    using energy3          = value<eV3_unit>;
    using inverse_energy   = value<inverse_eV_unit>;
    using inverse_energy2  = value<inverse_eV2_unit>;
    using inverse_energy3  = value<inverse_eV3_unit>;

    // set up some default units
    // note only long double allowed on a user-defined literal operator

    // eV is the fundamental unit
    constexpr value<eV_unit> operator "" _eV(long double d)
      {
        return energy(d);
      }

    // reduced Planck mass, expressed in eV
    constexpr value<eV_unit> operator "" _Mp(long double d)
      {
        return energy(2.436E27*d);
      }


    // overload arithmetic operators to allow dimensionful to be combined, compared

    template <typename Unit>
    constexpr value<Unit> operator+(const value<Unit>& a, const value<Unit>& b)
      {
        return value<Unit>(a.val + b.val);
      }

    template <typename Unit>
    constexpr bool operator<(const value<Unit>& a, const value<Unit>& b)
      {
        return(a.val < b.val);
      }

    template <typename Unit>
    constexpr bool operator>(const value<Unit>& a, const value<Unit>& b)
      {
        return(a.val > b.val);
      }


    // DIMENSIONFUL MULTIPLICATION

    constexpr energy2 operator*(const energy& a, const energy& b)
      {
        return energy2(a.val * b.val);
      }

    constexpr energy3 operator*(const energy& a, const energy2& b)
      {
        return energy3(a.val * b.val);
      }

    constexpr inverse_energy2 operator*(const inverse_energy& a, const inverse_energy& b)
      {
        return inverse_energy2(a.val * b.val);
      }

    constexpr inverse_energy3 operator*(const inverse_energy& a, const inverse_energy2& b)
      {
        return inverse_energy3(a.val * b.val);
      }

    constexpr inverse_energy3 operator*(const inverse_energy2& a, const inverse_energy& b)
      {
        return inverse_energy3(a.val * b.val);
      }


    // DIMENSIONFUL DIVISION

    constexpr energy operator/(const inverse_energy& a, const inverse_energy2& b)
      {
        return energy(a.val/b.val);
      }


    // DIMENSIONLESS RATIOS

    constexpr double operator/(const energy& a, const energy& b)
      {
        return a.val/b.val;
      }

    constexpr double operator/(const inverse_energy& a, const inverse_energy& b)
      {
        return a.val/b.val;
      }

    constexpr double operator*(const energy& a, const inverse_energy& b)
      {
        return a.val*b.val;
      }


    // DIMENSIONLESS MULTIPLICATION

    constexpr energy operator*(double a, const energy& b)
      {
        return energy(a * b.val);
      }

    constexpr energy2 operator*(double a, const energy2& b)
      {
        return energy2(a * b.val);
      }

    constexpr energy3 operator*(double a, const energy3& b)
      {
        return energy3(a * b.val);
      }

    constexpr inverse_energy operator*(double a, const inverse_energy& b)
      {
        return inverse_energy(a * b.val);
      }

    constexpr inverse_energy2 operator*(double a, const inverse_energy2& b)
      {
        return inverse_energy2(a * b.val);
      }

    constexpr inverse_energy3 operator*(double a, const inverse_energy3& b)
      {
        return inverse_energy3(a * b.val);
      }


    // DIMENSIONLESS DIVISION

    constexpr energy operator/(const energy& a, double b)
      {
        return energy(a.val / b);
      }

    constexpr energy2 operator/(const energy2& a, double b)
      {
        return energy2(a.val / b);
      }

    constexpr energy3 operator/(const energy3& a, double b)
      {
        return energy3(a.val / b);
      }

    constexpr inverse_energy operator/(const inverse_energy& a, double b)
      {
        return inverse_energy(a.val / b);
      }

    constexpr inverse_energy2 operator/(const inverse_energy2& a, double b)
      {
        return inverse_energy2(a.val / b);
      }

    constexpr inverse_energy3 operator/(const inverse_energy3& a, double b)
      {
        return inverse_energy3(a.val / b);
      }


    // RECIPROCALS

    constexpr energy operator/(double a, const inverse_energy& b)
      {
        return energy(a / b.val);
      }

    constexpr energy2 operator/(double a, const inverse_energy2& b)
      {
        return energy2(a / b.val);
      }

    constexpr inverse_energy operator/(double a, const energy& b)
      {
        return inverse_energy(a / b.val);
      }

    constexpr inverse_energy2 operator/(double a, const energy2& b)
      {
        return inverse_energy2(a / b.val);
      }


    // express SI-type units

    constexpr energy PlanckMass = 1.0_Mp;

    // numerical constant here is sqrt(1/8pi); note we have to write the literal explicitly
    // in C++11 (and probably C++14) because there is no constexpr square root function
    constexpr inverse_energy  sqrt_NewtonG = 0.1994711402007163 / PlanckMass;

    constexpr inverse_energy  Metre    = sqrt_NewtonG / 1.616199E-35;
    constexpr energy          Kilogram = 1.0 / (2.17651E-8 * sqrt_NewtonG);
    constexpr inverse_energy  Second   = sqrt_NewtonG / 5.39106E-44;
    constexpr energy          Kelvin   = 1.0 / (1.416833E32 * sqrt_NewtonG);

    constexpr inverse_energy  Kilometre = 1000 * Metre;
    constexpr inverse_energy  Mpc       = 3.08567758E22 * Metre;

    constexpr inverse_energy3 Mpc3     = Mpc*Mpc*Mpc;

    constexpr double          c        = 299792458 * Metre / Second;

  }



#endif //LSSEFT_UNITS_H
