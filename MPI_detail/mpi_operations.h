//
// Created by David Seery on 10/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_MPI_OPERATIONS_H
#define LSSEFT_MPI_OPERATIONS_H


namespace MPI_detail
  {

    constexpr unsigned int RANK_MASTER = 0;

    // message tags
    constexpr unsigned int MESSAGE_TERMINATE = 999;

  }   // namespace MPI


#endif //LSSEFT_MPI_OPERATIONS_H
