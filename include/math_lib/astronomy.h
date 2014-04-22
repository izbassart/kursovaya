#ifndef H_ASTRONOMY_H
#define H_ASTRONOMY_H

#include <algorithm>
#include <sstream>

#include "math_lib.h"
#include "quat.h"

namespace Astronomy
{

bool is_leap(int year);

int days_in_month(int month, int year);

int day_of_year(int day, int month, int year);

std :: pair<int, int> get_day_month(int doy, int year);


struct DateTime
{
  int year;
  int month;
  int date;
  
  int hours;
  int minutes;
  double seconds;


public:

  DateTime() : year(0), month(0), date(0), hours(0), minutes(0), seconds(0) { }; 
  DateTime(int y, int m, int d) : year(y), month(m), date(d), hours(0), minutes(0), seconds(0) { }
  DateTime(int y, int m, int d, int h, int mt, double s) : year(y), month(m), date(d), hours(h), minutes(mt), seconds(s) { }
  
  std :: string get_string() const;
  std :: string get_day_string() const;
  std :: string get_time_string() const;
 
  int DaysInThisYear() const
  {
    return is_leap(year) ? 366 : 365;
  }

  void AddDays(int d);
  void AddSeconds(double sec); // добавление секунд 
  bool operator<(const DateTime& dt) const;
  bool operator>(const DateTime& dt) const;
  bool operator==(const DateTime& dt) const;
  bool operator!=(const DateTime& dt) const
  {
    return !(*this == dt);
  }

  double operator-(const DateTime& dt) const;
  DateTime operator+(double T) const;
  DateTime operator-(double T) const;
  DateTime day_before() const;
  DateTime day_after() const;
  DateTime days_before(size_t d) const;
  DateTime days_after(size_t d) const;
  double TodaySeconds() const { return (hours * 60 + minutes) * 60 + seconds; }

  // задание даты в формате sql запроса oracle, генерирует TO_DATE
  std :: string oracle_sql_rep_date() const;
  std :: string oracle_sql_rep_datetime() const;
};

class JulianDate
{

public:
    long julian_days;

    JulianDate(long y, long m, long d)
    {
        julian_days = d - 32075 + 1461 * (y + 4800 + (m - 14) / 12) / 4 +
        367 * (m - 2 - ((m - 14) / 12) * 12) / 12 - 3 * ((y + 4900 + (m - 14) / 12) / 100) / 4;
    }

    JulianDate(long julian_days)
    {
        julian_days = julian_days;
    }

    double GetJD_days(int h, int m, double sec)
    {
        double th_s = 0.001 * sec + 0.06 * m + 3.6 * h;
        return julian_days - 0.5 /*5 / 8*/ + th_s / 86.4;
    }

    // time is returned in 1000 seconds
    double GetJD_sec(int h, int m, double sec)
    {
        double th_s = 0.001 * sec + 0.06 * m + 3.6 * h;
        return 86.4 * (julian_days - 0.5) + th_s;
    }

    long GetDays()
    {
        return julian_days;
    }

};

// время указывается в UTC
double DateTimeToJs(const DateTime& dt);

// время указывается в UTC
double DateTimeToJd(const DateTime& dt);

// время возвращается в UTC
DateTime JsToDateTime(double julianDate);

// время возвращается в UTC
inline DateTime JdToDateTime(double julianDate)
{
  return JsToDateTime(julianDate * 86.4);
}

// вычисление звездного времени
// время указывается в UTC
double GetStarTime(const DateTime& time);

// перевод из ИСК в ГСК
SVec3f ICStoGCS(const SVec3f& v, double starTime);

// перевод из ГСК в ИСК
SVec3f GCStoICS(const SVec3f& v, double starTime);

class SunPosition
{
  long julian_days;

public:
  SunPosition(int y, int m, int d)
  {
      JulianDate jd(y, m, d);
      julian_days = jd.GetDays();
  }

    /// Gets the sun position in ICS. See the AstronomyMath class to convert to GCS.
  
  SVec3f ComputeInICS(int h, int m, double s) const;
};


inline SVec3f GetSunPositionGreenwich(const DateTime& datetime)
{
  SVec3f p = SunPosition(datetime.year, datetime.month, datetime.date)
    .ComputeInICS(datetime.hours, datetime.minutes, datetime.seconds);
            
  return ICStoGCS(p, GetStarTime(datetime));
}

inline SVec3f GetSunPositionICS(const DateTime& datetime)
{
  return SunPosition(datetime.year, datetime.month, datetime.date)
    .ComputeInICS(datetime.hours, datetime.minutes, datetime.seconds);
}

// вычисляет кватернион рассогласования ОСК и ГСК
SQuat compute_orbital_SK(const SVec3f& position, const SVec3f& velocity);

// вычисляет кватерион рассогласования ИСК и ГСК
SQuat compute_inertial_SK(const Astronomy::DateTime& dt);

};



#endif