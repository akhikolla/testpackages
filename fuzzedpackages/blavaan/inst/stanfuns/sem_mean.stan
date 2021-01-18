  vector[] sem_mean(vector[] alpha, real[,,] B, real[,,] gamma, int[] g, int k, int Ng, int gamind, real[,] meanx){
    matrix[k,k] iden;
    vector[k] evlv[Ng];

    iden = diag_matrix(rep_vector(1.0, k));

    for(j in 1:Ng){
      if(gamind == 1){
        evlv[j] = inverse(iden - to_matrix(B[,,j])) * (alpha[j] + to_matrix(gamma[,,j]) * to_vector(meanx[,j]));

      } else {
        evlv[j] = inverse(iden - to_matrix(B[,,j])) * alpha[j];
      }
    }

    return evlv;
  }
