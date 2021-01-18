//

#include "lmbr.h"


void Clmbr::set_y( void )
// during initialization, set y values
{
	int i;
	for(i=0;i<n;i++) 
		if ( !R_FINITE( *(y_in+i) )  )  stop( _("invalid y value") );

	Vector<double>  y(n,0.),  irSy(n,0.);

	for (i=0;i<n;i++)  if( model_in > 0 )  y[i] = y_in[i];  else  y[i]= y_in[n-1-i];

	irSy = y;
	if( vectorS )  for(i=0;i<n;i++)  irSy[i] =  *( irS + i ) * y[i];
	if( matrixS )  for(i=0;i<n;i++) {
		irSy[i] = 0.;
		for( int k=0; k<n; k++ )  irSy[i] += *( irS + k*n + i ) * y[k];
	} 

// 'set_sy'  takes sy-values in original order
//  i.e. not reversed for 'model_in' = -2 or -3
	double*  irsy= Calloc( n, double );
	for (i=0;i<n;i++)
		if( model_in > 0 )  irsy[i] = irSy[i];  else  irsy[i]= irSy[n-1-i];

	set_sy( irsy, INIT );

	Free( irsy );

	return;
}

