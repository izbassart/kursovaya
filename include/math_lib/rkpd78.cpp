#include "rkdp78.h"

using namespace RKDP;

double  Solver78 ::A[14][13] = 
                       { 
                         { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
                         { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
                         { 0, 1.0 / 18.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
                         { 0, 1.0 / 48.0, 1.0 / 16.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
                         { 0, 1.0 / 32.0, 0, 3.0 / 32.0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
                         { 0, 5.0 / 16.0, 0, -75.0 / 64.0, 75.0 / 64.0, 0, 0, 0, 0, 0, 0, 0, 0 },
                         { 0, 3.0 / 80.0, 0, 0, 3.0 / 16.0, 3.0 / 20.0, 0, 0, 0, 0, 0, 0, 0 },
                         { 0, 2944.3841 / 61456.3906, 0, 0, 7773.6538 / 69253.8347, -2869.3883 / 112500.0, 2312.4283 / 180000.0, 
                           0, 0, 0, 0, 0, 0 },
                         { 0, 1601.6141 / 94669.2911, 0, 0, 6156.418 / 15873.2637, 2278.9713 / 63344.5777, 54581.5736 / 277105.7229,
                           -18019.3667 / 104330.7555, 0, 0, 0, 0, 0 },
                         { 0, 3963.2708 / 57359.1083, 0, 0, -43363.6366 / 68370.1615, -42173.9975 / 261629.2301, 10030.2831 / 72342.3059,
                           79020.4164 / 83981.3087, 80063.531 / 378307.1287, 0, 0, 0, 0 },
                         { 0, 24612.1993 / 134084.7787, 0, 0, -3769504.2795 / 1526876.6246, -30912.1744 / 106122.7803, 
                           -1299.2083 / 49076.6935, 600594.3493 / 210894.7869, 39300.6217 / 139667.3457, 12387.2331 / 100102.9789,
                           0, 0, 0 },
                         { 0, -102846.8189 / 84618.0014, 0, 0, 847823.5783 / 50851.2852, 131172.9495 / 143242.2823, -1030412.9995 / 170130.4382,
                           -4877792.5059 / 304793.956, 1533672.6248 / 103282.4649, -4544286.8181 / 339846.7696, 306599.3473 / 59717.2653, 0, 0 },
                         { 0, 18589.2177 / 71811.6043, 0, 0, -318509.4517 / 66710.7341, -47775.5414 / 109805.3517, -70363.5378 / 23073.9211,
                           573156.6787 / 102754.5527, 523286.6602 / 85006.6563, -409366.4535 / 80868.8257, 396213.7247 / 180595.7418,  6568.6358 / 48791.0083, 0 },
                         { 0, 40386.3854 / 49106.3109, 0, 0, -506849.2393 / 43474.0067, -41142.1997 / 54304.3805, 65278.3627 / 91429.6604,
                           1117396.2825 / 92532.0556, -1315899.0841 / 618472.7034, 393664.7629 / 197804.968, -16052.8059 / 68517.8525,
                           24863.8103 / 141353.1060, 0 }
                        };

double Solver78 :: B[] = {  0.0, // 0
                          1400.5451 / 33548.0064, // 1 
                          0.0, 0.0, 0.0, 0.0,  // 2, 3, 4, 5
                          -5923.8493 / 106827.7825, // 6
                          18160.6767 / 75886.7731, // 7
                          56129.2985 / 79784.5732, // 8 
                          -104189.143 / 137134.3529, // 9 
                          76041.7239 / 115116.5299, // 10 
                          11882.0643 / 75113.8087, // 11 
                          -52874.7749 / 222060.717, // 12
                          1.0 / 4  // 13
                        };

double Solver78 :: C[] =  { 
                          0.0, 0.0, // 0, 1
                          1.0 / 18, // 2
                          1.0 / 12, // 3
                          1.0 / 8, // 4
                          5.0 / 16, // 5
                          3.0 / 8, // 6
                          59.0 / 400, // 7
                          93.0 / 200, // 8
                          549002.3248 / 971916.9821, // 9 
                          13.0 / 20, // 10
                          120114.6811 / 129901.9798, // 11 
                          1.0, // 12
                          1.0  // 13
                        };

void Solver78 :: InitD()
{
  D[0] = 0;
  D[1] = 1345.1932 / 45517.6623;
  D[2] = D[3] = D[4] = D[5] = 0;
  D[6] = -80871.9846 / 97600.0145;
  D[7] = 175700.4468 / 564515.9321;
  D[8] = 65604.5339 / 26589.1186;
  D[9] = -386757.4721 / 151851.7206;
  D[10] = 46588.5868 / 32273.6535;
  D[11] = 5301.1238 / 66751.6719;
  D[12] = 2.0 / 45;
  D[13] = 0;
}

void Solver78 :: InitVar()
{
  InitD();
  for(int i = 1; i < 14; i++)
        D[i] = B[i] - D[i];

}

std :: vector<double>  Solver78 :: check_constants() 
{
  std :: vector<double> ret_value;
  
  InitD();

  double S;
  for(int i = 2; i < 14; i++)
  {
    S = -C[i];
    for (int j = 1; j < i; j++)
      S = S + A[i][j];
    ret_value.push_back(S);
  }

      
  for( int i = 3; i < 14; i++)
  {
    S = -(C[i] * C[i]) / 2;
    for (int j = 2 ;  j < i; j++)
      S = S + A[i][j] * C[j];
    ret_value.push_back(S);
  }
             
  for (int i=3; i < 14; i++)
  {
    S = C[i];
    S = - S * S * S / 3;
    for (int j = 2; j < i; j++) 
      S += A[i][j]*(C[j] * C[j]);
    ret_value.push_back(S);
  }
      
      
  for(int i = 6; i < 14; i++) 
  {
    S = C[i];
    S = -S * S * S * S / 4;
    for(int j = 2;  j < i; j++) 
      S += A[i][j] * C[j] * C[j] * C[j];
    ret_value.push_back(S);
  }
      
      
  for(int i = 1; i < 9; i++) 
  {
    S = - 1.0 / (i + 1);
    for(int j = 2; j < 14; j++)
      S += B[j] * exp(i * log(C[j]));
    ret_value.push_back(S);
  }
        
       

  for(int i = 1; i < 9; i++) 
  {
    S = -1.0 / (i + 1);
    for(int j = 2; j < 14; j++)
      S += D[j] * exp(i * log(C[j]));
    ret_value.push_back(S);
  }

      
  S = -1; 
  for(int i = 1; i < 14; i++) 
    S += B[i];
  ret_value.push_back(S);
     
  S = -1; 
  for(int i = 1; i < 14; i++)
    S += D[i];
  ret_value.push_back(S);
    
  InitVar();

  return ret_value;
}

bool Solver78 :: integrate(size_t Nvar, 
                         double * Y, 
                         double t_beg, 
                         double t_end, 
                         double beg_step, 
                         double start_precision,
                         bool control,
                         RightPart * in_rp,
                         OutPut * in_op)
{
  if (t_beg == t_end || beg_step <= 0 || start_precision <= 0|| control)
        return false;
      std :: vector<std :: vector<double> > F;
      F.resize(15);
      for(int q = 1 ; q < 15; q++)
        F[q].resize(Nvar);

      
      double TE;
          
            
      double T = t_beg;
      double step = T < t_end ? beg_step : -beg_step;
      while (fabs(t_end - T) >= 1e-11)
      {
        if( t_end > t_beg && (step > t_end - T) || t_end < t_beg && (step < t_end - T))
          step = t_end - T;
        in_op->DoOutput(T, step, Nvar, Y);
        if(control)
        {
          t_end = T;
          return true;
        }

        in_rp->Compute(Nvar, T, Y, &(F[1][0]));
        bool BOL = false;
        do
        {
          for(size_t i = 2; i < 14; i++)
          {
      
            for(size_t k = 0;  k < Nvar; k++)
            {
              double S = 0;
              for (size_t j = 1; j < i; j++)
                S += A[i][j] * F[j][k];
              F[14][k] = Y[k] + step * S;
            }
            in_rp->Compute(Nvar, T + step * C[i], &(F[14][0]), &(F[i][0]));
          }
          TE = 0;
          for(size_t  k = 0; k < Nvar; k++)
          {
            double TE1 = 0;
            for (int j = 1; j < 14; j++)
              TE1 = TE1 + D[j]*F[j][k];
            TE1 = fabs(TE1*step); 
            if (TE1 > TE) 
              TE = TE1;
          }

          if(TE > start_precision) 
          { 
            BOL = false; 
            step *= 0.7;
          }
          else 
            BOL = true;
        } 
        while(!BOL);

        for(size_t k = 0; k < Nvar; k++)
        {
          double S = 0;
          for(int i = 1; i < 14; i++)
            S += B[i] * F[i][k];
          Y[k] += S * step;
        }
    
        T += step;
        step *= 0.9 * exp(0.125 * log(start_precision / TE));
      }

      in_op->DoOutput(T, step, Nvar, Y); 
      return true;
}