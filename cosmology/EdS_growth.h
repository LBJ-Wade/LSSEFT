//
// Created by David Seery on 11/04/2017.
// --@@ // Copyright (c) 2017 University of Sussex. All rights reserved.
//
// This file is part of the Sussex Effective Field Theory for
// Large-Scale Structure platform (LSSEFT).
//
// LSSEFT is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// LSSEFT is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with LSSEFT.  If not, see <http://www.gnu.org/licenses/>.
//
// @license: GPL-2
// @contributor: David Seery <D.Seery@sussex.ac.uk>
// --@@
//

#ifndef LSSEFT_EDS_GROWTH_H
#define LSSEFT_EDS_GROWTH_H


// helper class to compute Einstein-de Sitter approximations to the growth
// functions and growth factors
template <typename ValueType>
class EdS_growth
  {
    
    // CONSTRUCTOR, DESTRUCTOR
    
  public:
    
    //! constructor captures linear growth factor and growth rate
    EdS_growth(ValueType D, ValueType f)
      : D_lin(D),
        D_lin_sq(D*D),
        D_lin_cb(D*D*D),
        f_lin(f)
      {
      }
    
    //! destructor is default
    ~EdS_growth() = default;
    
    
    // INTERFACE
    
  public:
    
    //! compute A
    ValueType DA() const { return 3.0*D_lin_sq/7.0; }
    
    //! compute B
    ValueType DB() const { return 2.0*D_lin_sq/7.0; }
    
    //! compute D
    ValueType DD() const { return 2.0*D_lin_cb/21.0; }
    
    //! compute E
    ValueType DE() const { return 4.0*D_lin_cb/63.0; }
    
    //! compute F
    ValueType DF() const { return 1.0*D_lin_cb/14.0; }
    
    //! compute G
    ValueType DG() const { return 1.0*D_lin_cb/21.0; }
    
    //! compute J
    ValueType DJ() const { return 1.0*D_lin_cb/9.0; }
    
    //! compute fA
    ValueType fA() const { return 2.0*f_lin; }
    
    //! compute fB
    ValueType fB() const { return 2.0*f_lin; }
    
    //! compute fD
    ValueType fD() const { return 3.0*f_lin; }
    
    //! compute fE
    ValueType fE() const { return 3.0*f_lin; }
    
    //! compute fF
    ValueType fF() const { return 3.0*f_lin; }
    
    //! compute fG
    ValueType fG() const { return 3.0*f_lin; }
    
    //! compute fJ
    ValueType fJ() const { return 3.0*f_lin; }
    
    
    // INTERNAL DATA
    
  private:
    
    //! capture linear growth function
    const ValueType D_lin;
    
    //! capture linear growth factor
    const ValueType f_lin;
    
    //! square linear growth function
    const ValueType D_lin_sq;
    
    //! cube linear growth function
    const ValueType D_lin_cb;
  
  };


#endif //LSSEFT_EDS_GROWTH_H
