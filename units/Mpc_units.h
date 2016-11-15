//
// Created by David Seery on 11/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_UNITS_H
#define LSSEFT_UNITS_H


#include "boost/math/constants/constants.hpp"
#include "boost/serialization/serialization.hpp"


namespace Mpc_units
  {

    // set up a template class representing a unit in a natural system where c=hbar=1
    // Here, we use Mpc as the base unit because this gives nice order-unity
    // values for most of our numbers

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

    using Mpc_unit          = unit<1>;
    using Mpc2_unit         = unit<2>;
    using Mpc3_unit         = unit<3>;
    using Mpc4_unit         = unit<4>;
    using inverse_Mpc_unit  = unit<-1>;
    using inverse_Mpc2_unit = unit<-2>;
    using inverse_Mpc3_unit = unit<-3>;
    using inverse_Mpc4_unit = unit<-4>;

    using energy           = value<inverse_Mpc_unit>;
    using energy2          = value<inverse_Mpc2_unit>;
    using energy3          = value<inverse_Mpc3_unit>;
    using energy4          = value<inverse_Mpc4_unit>;
    using inverse_energy   = value<Mpc_unit>;
    using inverse_energy2  = value<Mpc2_unit>;
    using inverse_energy3  = value<Mpc3_unit>;
    using inverse_energy4  = value<Mpc4_unit>;


    // set up some default units
    // note only long double allowed on a user-defined literal operator

    // Mpc is the fundamental unit
    constexpr value<Mpc_unit> operator "" _Mpc(long double d)
      {
        return inverse_energy(d);
      }


    // overload arithmetic operators to allow dimensionful to be combined, compared
    

    template <typename Unit>
    constexpr value<Unit> operator+(const value<Unit>& a, const value<Unit>& b)
      {
        return value<Unit>(a.val + b.val);
      }

    template <typename Unit>
    constexpr value<Unit> operator-(const value<Unit>& a, const value<Unit>& b)
      {
        return value<Unit>(a.val - b.val);
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

    constexpr energy3 operator*(const energy2& a, const energy& b)
      {
        return energy3(a.val * b.val);
      }

    constexpr energy4 operator*(const energy& a, const energy3& b)
      {
        return energy4(a.val * b.val);
      }

    constexpr energy4 operator*(const energy3& a, const energy& b)
      {
        return energy4(a.val * b.val);
      }

    constexpr energy4 operator*(const energy2& a, const energy2& b)
      {
        return energy4(a.val * b.val);
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

    constexpr inverse_energy2 operator*(const energy& a, const inverse_energy3& b)
      {
        return inverse_energy2(a.val * b.val);
      }

    constexpr inverse_energy2 operator*(const inverse_energy3& a, const energy& b)
      {
        return inverse_energy2(a.val * b.val);
      }

    constexpr inverse_energy operator*(const energy2& a, const inverse_energy3& b)
      {
        return inverse_energy(a.val * b.val);
      }

    constexpr inverse_energy operator*(const inverse_energy3& a, const energy2& b)
      {
        return inverse_energy(a.val * b.val);
      }

    constexpr inverse_energy4 operator*(const inverse_energy& a, const inverse_energy3& b)
      {
        return inverse_energy4(a.val * b.val);
      }

    constexpr inverse_energy4 operator*(const inverse_energy3& a, const inverse_energy& b)
      {
        return inverse_energy4(a.val * b.val);
      }

    constexpr inverse_energy3 operator*(const energy& a, const inverse_energy4& b)
      {
        return inverse_energy3(a.val * b.val);
      }

    constexpr inverse_energy3 operator*(const inverse_energy4& a, const energy& b)
      {
        return inverse_energy3(a.val * b.val);
      }


    // DIMENSIONFUL DIVISION

    constexpr energy operator/(const inverse_energy& a, const inverse_energy2& b)
      {
        return energy(a.val/b.val);
      }


    // DIMENSIONLESS RATIOS

    template <typename Unit>
    constexpr double operator/(const value<Unit>& a, const value<Unit>& b)
      {
        return a.val/b.val;
      }

    constexpr double operator*(const energy& a, const inverse_energy& b)
      {
        return a.val*b.val;
      }

    constexpr double operator*(const energy2& a, const inverse_energy2& b)
      {
        return a.val*b.val;
      }

    constexpr double operator*(const energy3& a, const inverse_energy3& b)
      {
        return a.val*b.val;
      }

    constexpr double operator*(const energy4& a, const inverse_energy4& b)
      {
        return a.val*b.val;
      }


    // DIMENSIONLESS MULTIPLICATION

    template <typename Unit>
    constexpr value<Unit> operator*(double a, const value<Unit>& b)
      {
        return value<Unit>(a * b.val);
      }

    template <typename Unit>
    constexpr value<Unit> operator*(const value<Unit>& a, double b)
      {
        return value<Unit>(a.val * b);
      }


    // DIMENSIONLESS DIVISION

    template <typename Unit>
    constexpr value<Unit> operator/(const value<Unit>& a, double b)
      {
        return value<Unit>(a.val / b);
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

    constexpr inverse_energy  Mpc          = 1.0_Mpc;
    constexpr inverse_energy2 Mpc2         = Mpc*Mpc;
    constexpr inverse_energy3 Mpc3         = Mpc*Mpc*Mpc;
    constexpr inverse_energy4 Mpc4         = Mpc*Mpc*Mpc*Mpc;

    constexpr inverse_energy  Metre        = Mpc / 3.08567758E22;
    constexpr inverse_energy  Kilometre    = 1000 * Metre;

    constexpr inverse_energy  sqrt_NewtonG = 1.616199E-35 * Metre;

    constexpr energy          Kilogram     = 1.0 / (2.17651E-8 * sqrt_NewtonG);
    constexpr inverse_energy  Second       = sqrt_NewtonG / 5.39106E-44;
    constexpr energy          Kelvin       = 1.0 / (1.416833E32 * sqrt_NewtonG);

    // numerical constant here is sqrt(1/8pi); note we have to write the literal explicitly
    // in C++11 (and probably C++14) because there is no constexpr square root function
    constexpr energy          PlanckMass   = 0.1994711402007163 / sqrt_NewtonG;
    constexpr energy          eV           = PlanckMass / 2.436E27;

    constexpr double          c            = 299792458 * Metre / Second;

  }   // namespace Mpc_units


// set up absolute value function
namespace std
  {
    
    template <typename Unit>
    constexpr Mpc_units::value<Unit> abs(const Mpc_units::value<Unit>& x)
      {
        return Mpc_units::value<Unit>(std::abs(x.val));
      }
    
  }   // namespace std



#endif //LSSEFT_UNITS_H
