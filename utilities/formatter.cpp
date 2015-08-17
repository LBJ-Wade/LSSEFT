//
// Created by David Seery on 13/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#include "formatter.h"


#include <iomanip>
#include <sstream>

#include "formatter.h"
#include "localizations/messages.h"


std::string format_time(boost::timer::nanosecond_type time, unsigned int precision)
	{
    std::ostringstream out;

    constexpr boost::timer::nanosecond_type mu_sec = 1000;
    constexpr boost::timer::nanosecond_type m_sec  = 1000*mu_sec;
    constexpr boost::timer::nanosecond_type sec    = 1000*m_sec;
    constexpr boost::timer::nanosecond_type minute = 60*sec;
    constexpr boost::timer::nanosecond_type hour   = 60*minute;

    if(time > hour)
	    {
        out << time/hour << FORMAT_HOUR_LABEL << " ";
        time = time % hour;
	    }
    if(time > minute)
	    {
        out << time/minute << FORMAT_MINUTE_LABEL << " ";
        time = time % minute;
	    }
    out << std::setprecision(precision) << static_cast<double>(time) / sec << FORMAT_SECOND_LABEL;

    return(out.str());
	}
