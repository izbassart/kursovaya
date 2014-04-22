#include "astronomy.h"
#include <iomanip>

using namespace Astronomy;

static int days_in_month_s[] =
{
  31, // January
  28, // February
  31, // March
  30, // April
  31, // May
  30, // June
  31, // July
  31, // August
  30, // September
  31, // October
  30, // November
  31  // December
};

static int days_by_month_usual[] =
{
  0, // January
  31, // February
  59, // March
  90, // April
  120, // May
  151, // June
  181, // July
  212, // August
  243, // September
  273, // October
  304, // November
  334  // December
};

static int days_by_month_leap[] =
{
  0, // January
  31, // February
  60, // March
  91, // April
  121, // May
  152, // June
  182, // July
  213, // August
  244, // September
  274, // October
  305, // November
  335  // December
};

const char* month_name[] = 
{
  "January",
  "February",
  "March",
  "April",
  "May",
  "June",
  "July",
  "August",
  "September",
  "October",
  "Novermber",
  "December"
};



bool Astronomy :: is_leap(int year) 
{
  return (year % 400 == 0) || ((year %  4 == 0) && (year % 100 != 0));
}

int days_in_month(int month, int year)
{
 if((month >= 1) && (month <= 12))
 {
   int days = days_in_month_s[month - 1];
   if((month == 2) && is_leap(year))
     days ++;
   return days;
 }
 return 0;
}

int Astronomy :: day_of_year(int day, int month, int year)
{
  int dd = day;
  if(is_leap(year))
  {
    dd += days_by_month_leap[month - 1];
  }
  else
  {
    dd += days_by_month_usual[month - 1];
  }
  return dd;
}

std :: pair<int, int> Astronomy :: get_day_month(int doy, int year)
{
  const int * dbm  = is_leap(year) ?  days_by_month_leap : days_by_month_usual;
  int m = 1;
  for( ; m < 12; m++)
    if(dbm[m] >= doy)
      break;
  return std :: pair<int, int>(doy - dbm[m - 1], m);
}

bool DateTime :: operator<(const DateTime& dt) const
{
  return DateTimeToJs(*this) < DateTimeToJs(dt);
}

bool DateTime :: operator>(const DateTime& dt) const
{
  return DateTimeToJs(*this) > DateTimeToJs(dt);
}

double DateTime :: operator-(const DateTime& dt) const
{
  return DateTimeToJs(*this) - DateTimeToJs(dt);
}

bool DateTime :: operator==(const DateTime& dt) const
{
  return year == dt.year && month == dt.month && date == dt.date && hours == dt.hours &&
         minutes == dt.minutes && seconds == dt.seconds;
}

DateTime DateTime :: operator+(double T) const
{
  return JsToDateTime(DateTimeToJs(*this) + T);
}

DateTime DateTime :: operator-(double T) const
{
  return JsToDateTime(DateTimeToJs(*this) - T);
}

void DateTime :: AddDays(int d)
{
  int cur_doy = day_of_year(date, month, year);
  cur_doy += d;
  int dity = DaysInThisYear();
  while(cur_doy > dity)
  {
    year ++;
    cur_doy -= dity;
    dity = DaysInThisYear();
  }

  std :: pair<int, int> dm = get_day_month(cur_doy, year);
  date = dm.first;
  month = dm.second;
}

void DateTime :: AddSeconds(double sec) // добавление секунд 
  {
    int d = int(sec / 86400);
    double sec_d = sec - d * 86400;
    int h = int(sec_d / 3600);
    double sec_dh = sec_d - h * 3600;
    int m = int(sec_dh / 60);
    double sec_dhm = sec_dh - m * 60;

    seconds += sec_dhm;
    if(seconds >= 60)
    {
      seconds -= 60;
      m++;
    }

    minutes += m;
    if(minutes >= 60)
    {
      minutes -= 60;
      h++;
    }

    hours += h;
    if(hours >= 24)
    {
      d++;
      hours -= 24;
    }

    AddDays(d);
  }

std :: string DateTime :: get_string() const
{
  std :: stringstream str;
  str << std ::setw(2) << date << "." << month_name[month - 1] << "." << year;
  str << " " << std ::setw(2) <<hours << ":" << std ::setw(2) << minutes << ":" << std ::setw(2) << seconds;
  return str.str();
}

std :: string DateTime :: oracle_sql_rep_date() const
{
  std :: stringstream str;
  str << "TO_DATE('" << date << "." << month << "." << year << "', 'DD.MM.YYYY')";
  return str.str();
}

std :: string DateTime :: oracle_sql_rep_datetime() const
{
  std :: stringstream str;
  str << "TO_DATE('" << date << "." << month << "." << year << " " << hours << ":" << minutes << ":" << int(seconds) << "', 'DD.MM.YYYY HH24:MI:SS')";
  return str.str();
}

std :: string DateTime :: get_day_string() const
{
  std :: stringstream str;
  str << std::setw(2) << date << "." << month_name[month - 1] << "." << year;
  return str.str();
}

std :: string DateTime :: get_time_string() const
{
  std :: stringstream str;
  str <<  " " << std::setw(2) << std :: setfill('0') << hours << ":" << 
         std ::setw(2) << std :: setfill('0') << minutes << ":" << std :: setfill('0') << std ::setw(2) << seconds;
  return str.str();
}
// возвращает момент времени на день раньше
DateTime DateTime :: day_before() const
{
  return *this - 86.4;
}


 DateTime DateTime :: days_before(size_t d) const
 {
   return *this - d * 86.4;
 }

 DateTime DateTime :: days_after(size_t d) const
 {
   return *this + d * 86.4;
 }

// возвращает момент времени на день позже
DateTime DateTime :: day_after() const
{
  return *this + 86.4;
}

double Astronomy :: DateTimeToJs(const DateTime& dt)
{
  JulianDate jd(dt.year, dt.month, dt.date);
  return jd.GetJD_sec(dt.hours, dt.minutes, dt.seconds);
}

double Astronomy :: DateTimeToJd(const DateTime& dt)
{
  JulianDate jd(dt.year, dt.month, dt.date);
  return jd.GetJD_days(dt.hours, dt.minutes, dt.seconds);
}

/*
DateTime Astronomy :: JsToDateTime(double julianDate)
{
  julianDate = julianDate / (86.4);
  DateTime date;
  double dblA, dblB, dblC, dblD, dblE, dblF;
  double dblZ, dblW, dblX;
  int day, month, year;
  dblZ = long(julianDate + 0.5);
  dblW = long((dblZ - 1867216.25) / 36524.25);
  dblX = long(dblW / 4);
  dblA = dblZ + 1 + dblW - dblX;
  dblB = dblA + 1524;
  dblC = long((dblB - 122.1) / 365.25);
  dblD = long(365.25 * dblC);
  dblE = long((dblB - dblD) / 30.6001);
  dblF = long(30.6001 * dblE);
  day = int(dblB - dblD - dblF);
  if (dblE > 13)
  {
      month = int(dblE - 13);
  }
  else
  {
      month = int(dblE - 1);
  }
  if ((month == 1) || (month == 2))
  {
      year = int(dblC - 4715);
  }
  else
  {
      year = int(dblC - 4716);
  }
  date = DateTime(year, month, day);

  JulianDate jd(year, month, day);

  date.AddSeconds(julianDate * (24 * 3.6) * 1000.0 - jd.GetDays() * 3600 * 24.0 + 43200);

  return date;
}
*/

DateTime Astronomy :: JsToDateTime(double julianDate)
{
  julianDate = julianDate / (86.4);
  DateTime date;
  
  long lZ = long(julianDate + 0.5);
  double dF = julianDate + 0.5 - lZ;
  long lA;
  if(lZ >= 2299161)
  {
    long lalpha = long((lZ - 1867216.25) / 36524.25);
    lA = lZ + 1 + lalpha - lalpha / 4;
  }
  else
    lA = lZ;

  long lB = lA + 1524;
  long lC = long((lB - 122.1) / 365.25);
  long lD = long(365.25 * lC);
  long lE = long((lB - lD) / 30.6001);

  double dd = lB - lD - long(30.6001 * lE) + dF;
  date.date = int(dd);

// Номер месяца есть:
  if(lE < 14)
  {
    date.month = lE - 1;
  }
  else
  {
    date.month = lE - 13;
  }

  if(date.month > 2)
  {
    date.year = lC - 4716;
  }
  else
  {
    date.year = lC - 4715;
  }

  date.AddSeconds(dF * 86400.0);

  return date;
}

double Astronomy :: GetStarTime(const DateTime& time)
{
  JulianDate jd(time.year, time.month, time.date);
  long julian_days = jd.GetDays();

  double s = time.seconds;
  double DJ = julian_days - 2415020.5;
  double TIN = 0.001 * s + 0.06 * time.minutes + 3.6 * time.hours;
  double SZ0 = 0.276919398 + 2.73790926493e-3 * DJ +
                1.075231E-6 * (DJ / 36525) * (DJ / 36525) +
                1.0027379093 * (TIN / 86.4 /*- 0.125*/);
  return 2 * M_PI * (SZ0 - (int)(SZ0));
}

SVec3f SunPosition :: ComputeInICS(int h, int m, double s) const
{
  double t = (60 * h + m) * 60 + s;
  double grad = M_PI / 180.0;
      
  double TE = (julian_days + t / 86400 - 2415020 - 0.5) / 36525;
  double angCE = grad * (23.452294 - 0.0130125 * TE);
  double SE = sin(angCE);
  double CE = cos(angCE);
  double EE = 0.01675104 - 0.0000418 * TE;
  double angLS = grad * (279.69668 + TE * (36000.76892 + TE * 0.0003025));
  double angMS = grad * (358.475833 + TE * (35999.04975 - TE * 0.000150));
  angLS += 2 * EE * sin(angMS) + 1.25 * EE * EE * sin(2 * angMS);
  double X = cos(angLS);
  double LS = sin(angLS);
  double Y = CE * LS;
  double Z = SE * LS;
  return SVec3f(X, Y, Z);
}

SVec3f Astronomy :: ICStoGCS(const SVec3f& v, double starTime)
{
  double cw = cos(starTime);
  double sw = sin(starTime);
  
  return SVec3f(
                cw * v.x + sw * v.y,
                -sw * v.x + cw * v.y,
                v.z);
}

SVec3f Astronomy :: GCStoICS(const SVec3f& v, double starTime)
{
  double cw = cos(starTime);
  double sw = sin(starTime);
  
  return SVec3f(cw * v.x + -sw * v.y,
                sw * v.x + cw * v.y,
                v.z);
}

// вычисляет кватернион рассогласования ОСК и ГСК
SQuat Astronomy :: compute_orbital_SK(const SVec3f& position, const SVec3f& velocity)
{
  
  // вычисление ортов ОСК OX, OY, OZ 

  SVec3f OY(position); 
  OY.Normalize();
  SVec3f OZ = velocity ^ OY;
  OZ.Normalize();
  SVec3f OX = OY ^ OZ;
  
  return gen_from_matrix(SMatrix3f(OX, OY, OZ));
}

// вычисляет кватерион рассогласования ИСК и ГСК
SQuat Astronomy :: compute_inertial_SK(const Astronomy::DateTime& dt)
{
  double sz2 = GetStarTime(dt) / 2;
  return SQuat(cos(sz2), 0, 0, sin(sz2));
}