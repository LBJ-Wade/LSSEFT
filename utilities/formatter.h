//
// Created by David Seery on 13/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_FORMATTER_H
#define LSSEFT_FORMATTER_H

#include <boost/timer/timer.hpp>


std::string format_time(boost::timer::nanosecond_type time, unsigned int precision=3);


#endif //LSSEFT_FORMATTER_H
