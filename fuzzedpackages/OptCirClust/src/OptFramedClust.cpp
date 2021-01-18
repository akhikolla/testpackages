//
//  OptFramedClust.cpp
//  Optimal framed clustering in polylog-linear time
//
//  Created by Tathagata Debnath on 6/15/20.
//  Copyright © 2020 Tathagata Debnath and Joe Song.
//  All rights reserved.
//  Revised by Joe Song.

#include "OptFramedClust.h"



clustering MFC(std::vector<double> & Data_Points,
               int width, int K,
               int First, int Last,
               int Prev, int Next)

{
  struct clustering cluster;

  if(Last >= First)
  {


    std::sort (Data_Points.begin(), Data_Points.end());

    int N = Data_Points.size();



    std::vector<std::vector<int> > Cluster_Border ( (Last  + 1), std::vector<int>(K,0) );

    std::vector<double> sum_x(N), sum_x_sq(N);



    double shift = Data_Points[N/2];

    double scale = std::max(abs(Data_Points[0]), abs(Data_Points[N-1]));


    if(scale == 0)
    {
      scale = 1;
    }

    sum_x[0] = (Data_Points[0] - shift)/scale;
    sum_x_sq[0] = (Data_Points[0] - shift) * (Data_Points[0] - shift)/( scale * scale );



    for(int i = 1; i < N; ++i) {

      sum_x[i] = sum_x[i-1] + (Data_Points[i] - shift)/scale;
      sum_x_sq[i] = sum_x_sq[i-1] + (Data_Points[i] - shift) * (Data_Points[i] - shift) / (scale * scale);
    }


    std::vector< std::vector< double > > S( K, std::vector<double>(width,0) );

    std::vector< std::vector< int > > J( K, std::vector<int>(width,0) );



    frame_info frame = BDP(
      width,  K, First,  Last, Prev,
      Next, S, J, sum_x, sum_x_sq, Cluster_Border);





    std::vector<double> withinss(K);






    cluster.Frame_ID = frame.Frame_ID;

    cluster.ssq = frame.ssq * scale * scale;

    cluster.Borders = Cluster_Border[frame.Frame_ID];

    cluster.centers = withinss;

    cluster.size = withinss;


    //   cluster.centers[0] = (cluster.Frame_ID + cluster.Borders[0])/2 ;

    cluster.size[0] = cluster.Borders[0] - cluster.Frame_ID + 1;

    if(cluster.Frame_ID > 0)
    {
      cluster.centers[0] = (scale * ( sum_x[cluster.Borders[0]] - sum_x[cluster.Frame_ID - 1])/cluster.size[0]) + shift;
    }else{
      cluster.centers[0] = (scale * sum_x[cluster.Borders[0]]  / cluster.size[0]) + shift ;
    }


    withinss[0] = ssq( cluster.Frame_ID, cluster.Borders[0], 0, sum_x, sum_x_sq )* scale * scale;

    for(int i=1; i < K; i++)
    {

      cluster.size[i] = cluster.Borders[i] - cluster.Borders[i-1];


      cluster.centers[i] = (scale * (sum_x[cluster.Borders[i]] - sum_x[cluster.Borders[i-1]]) /cluster.size[i])  + shift ;


      withinss[i] = ssq( cluster.Borders[i-1]+1, cluster.Borders[i], 0, sum_x, sum_x_sq )* scale * scale;

    }


    cluster.totss = ssq( cluster.Frame_ID, cluster.Borders[K-1], 0, sum_x, sum_x_sq ) * scale * scale;

    cluster.withinss = withinss;


    return cluster;

  }
  else
  {
    return cluster;
  }
}



frame_info BDP( int width, int K,
                int First, int Last,
                int Prev, int Next,
                std::vector< std::vector< double > > &  S,
                std::vector< std::vector< int > > & J,
                const std::vector<double> & sum_x,
                const std::vector<double> & sum_x_sq,
                std::vector< std::vector<int> > & Cluster_Border)
{
  struct frame_info frame_mid ;

  if(Last >= First )

  {

    int Middle_Frame = (First + Last) / 2;

    linear_clustering(S, J, Prev, Next, Middle_Frame, sum_x, sum_x_sq, Cluster_Border);



    frame_mid.ssq = S[K-1][width-1];

    frame_mid.Frame_ID = Middle_Frame;

    int j = width - 1;

    int k = K - 1;

    Cluster_Border[Middle_Frame][k] = j + Middle_Frame;

    while(k > 0)
    {

      j = J[k][j] - 1;

      k = k - 1;

      Cluster_Border[Middle_Frame][k] = j + Middle_Frame;



    }







    frame_info frame_left = BDP(
      width,  K, First,  Middle_Frame - 1, Prev,  Middle_Frame,
      S, J, sum_x, sum_x_sq, Cluster_Border
    );

    frame_info frame_right = BDP(
      width,  K,  Middle_Frame + 1, Last, Middle_Frame, Next,
      S, J,  sum_x, sum_x_sq, Cluster_Border
    );





    if(frame_mid.ssq < frame_left.ssq && frame_mid.ssq <= frame_right.ssq)
    {
      return  frame_mid;
    }
    else if(frame_left.ssq <= frame_mid.ssq && frame_left.ssq  <= frame_right.ssq)
    {
      return  frame_left;
    }
    else
    {
      return  frame_right;
    }
  }
  else{
    return frame_mid;
  }


}








////////////////////////////////////////////////////////////////////////////////////////////////////////////


void linear_clustering(std::vector< std::vector< double > > & S,
                       std::vector< std::vector< int > > & J,
                       int Prev, int Next, int Middle_Frame,
                       const std::vector<double> & sum_x,
                       const std::vector<double> & sum_x_sq,
                       std::vector< std::vector<int> > & Cluster_Border
)
  /*
   S: K x N matrix. S[q][i] is the sum of squares of the distance from
   each x[i] to its cluster mean when there are exactly x[i] is the
   last point in cluster q
   J: K x N backtrack matrix

   NOTE: All vector indices in this program start at position 0
   */
{

  const int K = (int) S.size();
  const int N = (int) S[0].size();



  int jmin, jmax, imin, imax;

  if(Prev > -1 )
  {
    imin =  std::max(0, (Cluster_Border[Prev][0] - Middle_Frame  ));
  }
  else
  {
    imin = 0;
  }

  if(Next > -1 )
  {
    imax =  std::min(N-1, (Cluster_Border[Next][0] - Middle_Frame  ));
  }
  else
  {
    imax= N-1;
  }

  for(int i = imin; i <= imax; ++i) {
    S[0][i] = ssq( 0, i, Middle_Frame, sum_x, sum_x_sq);
  }


  for(int k = 1; k < K; ++k) {


    if(k < K - 1) {

      if(Prev > -1 )
      {

        jmin = std::max(k, (Cluster_Border[Prev][k-1] - Middle_Frame + 1) );

        imin =  std::max(k, (Cluster_Border[Prev][k] - Middle_Frame  ));


      }
      else
      {
        jmin = k;

        imin = std::max(1, k);
      }

      if(Next > -1 )
      {
        jmax = std::min(N-1, (Cluster_Border[Next][k-1] - Middle_Frame + 1) );

        imax = std::min(N-1, (Cluster_Border[Next][k] - Middle_Frame  ));

      }
      else
      {
        jmax = N-1;

        imax= N-1;
      }




    } else {
      // No need to compute S[K-1][0] ... S[K-1][N-2]

      if(Prev > -1 )
      {
        jmin = std::max(k, (Cluster_Border[Prev][k-1] - Middle_Frame + 1) );
      }
      else
      {
        jmin = k;
      }

      if(Next > -1 )
      {
        jmax = std::min(N-1, (Cluster_Border[Next][k-1] - Middle_Frame + 1) );
      }
      else
      {
        jmax = N-1;
      }

      imax = imin = N-1;
    }




    fill_row_k(imin, imax, k, Middle_Frame, jmin, jmax, S, J, sum_x, sum_x_sq);



  }



}




/////////////////////////////////////////////////////////////////////////////////////////////////////////////


/*
 void fill_row_q_2020_07_22(int imin, int imax, int q, int imin_0, int imax_0, int Middle_Frame,
 int jmin, int jmax,
 std::vector< std::vector<double> > & S,
 std::vector< std::vector<int> > & J,
 const std::vector<double> & sum_x,
 const std::vector<double> & sum_x_sq)
 {
 if(imin > imax) {
 return;
 }

 const int N = (int) S[0].size();

 int i = (imin + imax) / 2;



 // Initialization of S[q][i]:

 // Check this portion

 // initialize  S[q][i] =  INF.

 S[q][i] = 1.79769e+308;
 J[q][i] = i;

 int jlow=q; // the lower end for j

 if(imin > q) {
 // jlow = std::max(jlow, (int)J[q][imin-1]);
 jlow = std::max(jlow, jmin);
 }
 // jlow = std::max(jlow, J[q-1][i]);

 int jhigh = i; // the upper end for j
 if(imax < N-1) {
 // jhigh = std::min(jhigh, (int)J[q][imax+1]);
 jhigh = std::min(jhigh, jmax);
 }

 for(int j=jhigh; j>=jlow; --j) {

 // compute s(j,i)
 double sji = ssq(j, i, Middle_Frame, sum_x, sum_x_sq);


 if(sji + S[q-1][jlow-1] >= S[q][i]) break;

 // Examine the lower bound of the cluster border
 // compute s(jlow, i)
 double sjlowi =
 ssq(jlow, i, Middle_Frame, sum_x, sum_x_sq);
 // ssq(jlow, i, sum_x, sum_x_sq, sum_w);

 double SSQ_jlow = sjlowi + S[q-1][jlow-1];

 if(SSQ_jlow < S[q][i]) {
 // shrink the lower bound
 S[q][i] = SSQ_jlow;
 J[q][i] = jlow;
 }
 jlow ++;

 double SSQ_j = sji + S[q - 1][j - 1];
 if(SSQ_j < S[q][i]) {
 S[q][i] = SSQ_j;
 J[q][i] = j;
 }
 }

 // int t;



 // jmin= (imin > q) ? J[q][imin-1] : q;

 int jmin_left = (imin > imin_0) ? J[q][imin-1] : jmin;

 int jmax_left = J[q][i];

 fill_row_q_2020_07_22(imin, i-1, q, imin_0, imax_0, Middle_Frame, jmin_left, jmax_left,
 S, J, sum_x, sum_x_sq);



 int jmin_right = J[q][i];

 // jmax = (imax < N-1) ? J[q][imax+1] : imax;

 int jmax_right = (imax < imax_0) ? J[q][imax+1] : jmax;

 fill_row_q_2020_07_22(i+1, imax, q, imin_0, imax_0, Middle_Frame, jmin_right, jmax_right,
 S, J, sum_x, sum_x_sq);






 }

 */









void fill_row_k(int imin, int imax, int k, int Middle_Frame,
                int jmin, int jmax,
                std::vector< std::vector<double> > & S,
                std::vector< std::vector<int> > & J,
                const std::vector<double> & sum_x,
                const std::vector<double> & sum_x_sq)
{
  if(imin > imax) {
    return;
  }

  // const int N = (int) S[0].size();

  int i = (imin + imax) / 2;



  // Initialization of S[q][i]:

  // Check this portion

  // initialize  S[q][i] =  INF.

  S[k][i] = std::numeric_limits<double>::infinity(); // 1.79769e+308;
  J[k][i] = i;

  int jlow=k; // the lower end for j

  if(imin > k) {
    // jlow = std::max(jlow, (int)J[q][imin-1]);
    jlow = std::max(jlow, jmin);
  }
  // jlow = std::max(jlow, J[q-1][i]);

  int jhigh = i; // the upper end for j
  /* if(imax < N-1) {
   // jhigh = std::min(jhigh, (int)J[q][imax+1]);
   jhigh = std::min(jhigh, jmax);
  } */

  jhigh = std::min(jhigh, jmax);

  for(int j=jhigh; j>=jlow; --j) {

    // compute s(j,i)
    double sji = ssq(j, i, Middle_Frame, sum_x, sum_x_sq);


    if(sji + S[k-1][jlow-1] >= S[k][i]) break;




    double SSQ_j = sji + S[k - 1][j - 1];
    if(SSQ_j < S[k][i]) {
      S[k][i] = SSQ_j;
      J[k][i] = j;
    }
  }



  fill_row_k(imin, i-1, k,  Middle_Frame, jmin, J[k][i],
             S, J, sum_x, sum_x_sq);





  fill_row_k(i+1, imax, k,  Middle_Frame, J[k][i], jmax,
             S, J, sum_x, sum_x_sq);






}



////////////////////////////////////////////////////////////////////////////////////////////////////////////



double ssq(const int j, const int i, const int Middle_Frame,
           const std::vector<double> & sum_x,
           const std::vector<double> & sum_x_sq
)
{
  double sji(0.0);


  if(j >= i) {
    sji = 0.0;
  } else if(j > 0 || Middle_Frame > 0) {
    double muji = (sum_x[i + Middle_Frame] - sum_x[j-1 + Middle_Frame]) / (i - j + 1);
    sji = sum_x_sq[i + Middle_Frame] - sum_x_sq[j-1 + Middle_Frame] - (i - j + 1) * muji * muji;



  } else {
    // j == 0 && Middle_Frame == 0

    sji = (sum_x_sq[i ] ) - (sum_x[i]) * (sum_x[i]) / (i+1);

  }


  sji = (sji < 0) ? 0 : sji;
  return sji;
}

////////////////////////////////////////////

void backtrack (std::vector< std::vector<int> > & J,
                std::vector<int> & B,
                int K, int N)
{
  int j = N - 1;

  int k = K - 1;

  B[k] = j;

  while(k > 0)
  {

    j = J[k][j] - 1;

    k = k - 1;

    B[k] = j;

  }

  return;

}



///////////////////////////////////////////////////////////////////////////////////////////////////




