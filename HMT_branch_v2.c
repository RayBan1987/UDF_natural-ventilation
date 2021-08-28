/*******************************************************************
UDF for specifying an energy source term to simulate fire source
*******************************************************************/
/*repair the bug of lf higher than 0.45m,20210816*/

#include "udf.h"

#define a 0.2
#define b 0.1
#define Q 10.04
#define surfNO 34 /*branch tunnel end surface ID, must be checked before calculating*/
#define Rg 287  /*universal gas constant, J/kgK*/

DEFINE_SOURCE(heat_source,c,t,dS,eqn)
{
    real xc[ND_ND];
    real D,lf,xmax,xmin,zmin,zmax,source,volume,Q_xing;
    C_CENTROID(xc,c,t);
    D=2*(a+b)/M_PI;   /*D is perimeter diamter*/
	Q_xing=Q/(1.2*1.0*300*sqrt(9.81)*pow(D,2.5));
	lf=(3.7*pow(Q_xing,0.4)-1.02)*D;  /*calculate the height of fire*/
	lf=(lf>0.45)?0.45:lf;
	
	volume=a*b*lf;/*calculate the volume of fire source*/
	
    xmax=0.3+a/2;
    xmin=0.3-a/2;
    zmax=4.5+b/2;
    zmin=4.5-b/2;
    
    /*now begins to judge*/
    if (xc[0]<=xmax && xc[0]>=xmin && xc[2]<=zmax && xc[2]>=zmin && xc[1]<=lf)
	{
		source=Q*1000/volume;
        dS[eqn] = 0.0;
	}
	else
	{
		source=0;
        dS[eqn] = 0.0;
	}
    return source;
}

DEFINE_SOURCE(mass_source,c,t,dS,eqn)
{
    real xc[ND_ND];
    real D,lf,xmax,xmin,zmin,zmax,source,volume,Q_xing;
    C_CENTROID(xc,c,t);
    D=2*(a+b)/M_PI;   /*D is perimeter diamter*/
	Q_xing=Q/(1.2*1.0*300*sqrt(9.81)*pow(D,2.5));
	lf=(3.7*pow(Q_xing,0.4)-1.02)*D;  /*calculate the height of fire*/
	lf=(lf>0.45)?0.45:lf;
	volume=a*b*lf;/*calculate the volume of fire source*/
	
    xmax=0.3+a/2;
    xmin=0.3-a/2;
    zmax=4.5+b/2;
    zmin=4.5-b/2;
    
    /*now begins to judge*/
    if (xc[0]<=xmax && xc[0]>=xmin && xc[2]<=zmax && xc[2]>=zmin && xc[1]<=lf)
	{
		if ((3.7*pow(Q_xing,0.4)-1.02)*D>0.45) /*flame touch the ceiling*/
		    {source=0.426*pow(Q,1.0/3)*pow(lf,2.0/3)/a/pow(b,3.0)*pow(xc[2]-4.5,2.0)+0.0355*pow(Q,1.0/3)*pow(lf,2.0/3)/(a*b);
		     dS[eqn] = 0.0;}
		else /*flame not touch the ceiling*/
			{source=0.1704*pow(Q,1.0/3)*pow(lf,2.0/3)/a/pow(b,3.0)*pow(xc[2]-4.5,2.0)+0.0568*pow(Q,1.0/3)*pow(lf,2.0/3)/(a*b);
		     dS[eqn] = 0.0;}
	}
	else
	{
		source=0;
        dS[eqn] = 0.0;
	}
    return source;
}


/*******************************************************************
UDF for recording the mass flux at portal Z=0
*******************************************************************/
DEFINE_INIT(head_z0,d)
{
	Thread *t, *tf;
	cell_t c;
	face_t f;
	real xc[ND_ND];
	real xf[ND_ND];
	real A[ND_ND];  /*store face normal vector*/
	int n=-1,cell_num=0; /*loop variable*/
	
	FILE *fp; /*pointer to a file */
	fp = fopen("flux_z0.txt","w+");
	fprintf(fp,"/*this is recorded by FLUENT UDF*/\n");
	
	/*get the face normal area and record*/
	thread_loop_c(t,d)
	{
		begin_c_loop(c,t)
		{
		  C_CENTROID(xc,c,t);  /*get the center coordinates of one cell*/
		  if (xc[2]<=0.5 && xc[2]>=-0.5) /*first, shrink search range*/
		  {
			   c_face_loop(c,t,n) /*loop over all the faces of one cell*/
			   {
			      f=C_FACE(c,t,n);
				  tf=C_FACE_THREAD(c,t,n);
				  F_CENTROID(xf,f,tf); /*get the coordinate of a face*/
				  
				 if ((fabs(xf[2]-0)<1e-4) && (xc[2]>0))/*indicates it is  boundary*/
				  {
				   /*here puts integeral*/
				   F_AREA(A,f,tf);
				   cell_num++;
				   /*fprintf(fp,"x=%g,y=%g,z=%g\n",xc[0],xc[1],xc[2]);
				   	fprintf(fp,"x=%g,y=%g,z=%g\n",xf[0],xf[1],xf[2]);*/
				  }
		       }
		   }
	    }end_c_loop(c,t)
	}
	fprintf(fp,"face normal area is %g, %g, %g,\t",A[0]/NV_MAG(A),A[1]/NV_MAG(A),A[2]/NV_MAG(A));
	fprintf(fp,"cell_num is %d\t",cell_num);
	if (ND_DOT(A[0]/NV_MAG(A),A[1]/NV_MAG(A),A[2]/NV_MAG(A),0,0,1)>0)
	{
		fprintf(fp,"face area normal points into the tunnel\n");
	}
	else
	{
		fprintf(fp,"face area normal points out of the tunnel\n");
	}
	fprintf(fp,"vf_out mf_out vf_in mf_in mass_flow_in_rate neutral_height q_in q_out cell_num\n");
	fclose(fp);
}

DEFINE_EXECUTE_AT_END(flux_z0)
{
	Domain *d;
   	Thread *t,*tf;
    cell_t c;
	face_t f;
	real xc[ND_ND];
	real xf[ND_ND];
	real A[ND_ND];  /*store face normal vector*/

	real d_vf=0, d_mf=0;
	real vf_in=0, mf_in=0;
	real vf_out=0, mf_out=0;
	real q_in=0, q_out=0;
	real cell_num=0;
	
	int n=-1; /*loop variable*/
	real neutral_height=0,A_orintation=0,epsilon=100; /*epsilon was used to judge pressure*/
	
	FILE *fp; /*pointer to a file */
	fp = fopen("flux_z0.txt","a+");
	d = Get_Domain(1);	
		
	thread_loop_c(t,d)
	{
		begin_c_loop(c,t)
		{
		  C_CENTROID(xc,c,t);  /*get the center coordinates of one cell*/
		  if (xc[2]<=0.5 && xc[2]>=0) /*first, shrink search range*/
		  {
			   c_face_loop(c,t,n) /*loop over all the faces of one cell*/
			   {
			      f=C_FACE(c,t,n);
				  tf=C_FACE_THREAD(c,t,n);
				  F_CENTROID(xf,f,tf); /*get the coordinate of a face*/
				  
				 if ((fabs(xf[2]-0)<1e-4) && (xc[2]>0))/*indicates it is  boundary*/
				  {
				   /*here puts integeral*/
				   F_AREA(A,f,tf);
				   /*  A_orintation=ND_DOT(A[0]/NV_MAG(A),A[1]/NV_MAG(A),A[2]/NV_MAG(A),0,0,1);/*to judge the orientation of A*/
				   d_vf=ND_DOT(A[0],A[1],A[2],C_U(c,t),C_V(c,t),C_W(c,t));/*dot product to calculate the volume flux*/
				    if (ND_DOT(C_U(c,t),C_V(c,t),C_W(c,t),0,0,1)<0)/*indicates area normal and vlocity orients the opposite, inflow*/
				    { /*outflow*/
					    vf_out=vf_out+fabs(d_vf);
						mf_out=mf_out+fabs(d_vf)*C_R(c,t);
						/*fprintf(fp,"out %-5g\n",C_W(c,t) );  write the coordinate*/
						q_out=q_out+C_CP(c,t)*C_R(c,t)*(C_T(c,t)-300)*fabs(d_vf);/*sum the heat flux*/
						 if ((fabs(xf[0]-0.3)<2e-1) && (epsilon <= fabs(C_P(c,t))))/*x= 0.1m-0.5m*/
					        {  
							   epsilon = fabs(C_P(c,t));
							   neutral_height=xf[1];
							}
							
					}
				  
					else 
						{
					     vf_in=vf_in+fabs(d_vf);
					     mf_in=mf_in+fabs(d_vf)*C_R(c,t);
					  /* fprintf(fp,"in %-5g\n",C_W(c,t) );  write the coordinate*/
					     q_in=q_in+C_CP(c,t)*C_R(c,t)*(C_T(c,t)-300)*fabs(d_vf);/*sum the heat flux*/
					   }
					  cell_num++;
				  }
		       }
		   }
	    }
		end_c_loop(c,t)  
	}
	fprintf(fp,"%-5g %-5g %-5g %-5g %-5g %-5g %-5g %-5g %-5g\n", vf_out,mf_out,vf_in,mf_in,mf_in-mf_out,neutral_height,q_in,q_out,cell_num);  /*write the coordinate*/
	fclose(fp); 
}

DEFINE_INIT(head_z9,d)
{
	Thread *t, *tf;
	cell_t c;
	face_t f;
	real xc[ND_ND];
	real xf[ND_ND];
	real A[ND_ND];  /*store face normal vector*/
	int n=-1,cell_num=0; /*loop variable*/
	
	FILE *fp; /*pointer to a file */
	fp = fopen("flux_z9.txt","w+");
	fprintf(fp,"/*this is recorded by FLUENT UDF*/\n");
	
	
	/*get the face normal area and record*/
	thread_loop_c(t,d)
	{
		begin_c_loop(c,t)
		{
		  C_CENTROID(xc,c,t);  /*get the center coordinates of one cell*/
		  if (xc[2]<=9 && xc[2]>=8.5) /*first, shrink search range*/
		  {
			   c_face_loop(c,t,n) /*loop over all the faces of one cell*/
			   {
			      f=C_FACE(c,t,n);
				  tf=C_FACE_THREAD(c,t,n);
				  F_CENTROID(xf,f,tf); /*get the coordinate of a face*/
				  
				 if ((fabs(xf[2]-9)<1e-3) && (xc[2]<=9)) /*indicates it is  boundary*/
				  {
				   /*here puts integeral*/
				   F_AREA(A,f,tf);
				   cell_num++;
				   /* fprintf(fp,"x=%g,y=%g,z=%g\n",xc[0],xc[1],xc[2]);
				   	fprintf(fp,"x=%g,y=%g,z=%g\n",xf[0],xf[1],xf[2]);*/
				  }
		       }
		   }
	    }end_c_loop(c,t)
	}
	fprintf(fp,"face normal area is %g, %g, %g, \t",A[0]/NV_MAG(A),A[1]/NV_MAG(A),A[2]/NV_MAG(A));
	fprintf(fp,"cell_num is %d\t",cell_num);
	if (ND_DOT(A[0]/NV_MAG(A),A[1]/NV_MAG(A),A[2]/NV_MAG(A),0,0,1)>0)
	{
			fprintf(fp,"face area normal points out the tunnel\n");
	}	
	else
	{
			fprintf(fp,"face area normal points in the tunnel\n");
	}
	fprintf(fp,"vf_out mf_out vf_in mf_in mass_flow_in_rate neutral_height q_in q_out cell_num\n");
	fclose(fp);
}

DEFINE_EXECUTE_AT_END(flux_z9)
{
	Domain *d;
   	Thread *t,*tf;
    cell_t c;
	face_t f;
	real xc[ND_ND];
	real xf[ND_ND];
	real A[ND_ND];  /*store face normal vector*/

	real d_vf=0, d_mf=0;
	real vf_in=0, mf_in=0;
	real vf_out=0, mf_out=0;
	real q_in=0, q_out=0;
	real cell_num=0;
	
	int n=-1; /*loop variable*/
	real neutral_height=0,A_orintation=0,epsilon=100;
	
	FILE *fp; /*pointer to a file */
	fp = fopen("flux_z9.txt","a+");
	d = Get_Domain(1);	
		
	thread_loop_c(t,d)
	{
		begin_c_loop(c,t)
		{
		  C_CENTROID(xc,c,t);  /*get the center coordinates of one cell*/
		  if (xc[2]<=9 && xc[2]>=8.5) /*first, shrink search range*/
		  {
			   c_face_loop(c,t,n) /*loop over all the faces of one cell*/
			   {
			      f=C_FACE(c,t,n);
				  tf=C_FACE_THREAD(c,t,n);
				  F_CENTROID(xf,f,tf); /*get the coordinate of a face*/
				  
				  if ((fabs(xf[2]-9)<1e-3) && (xc[2]<=9)) /*indicates it is  boundary*/
				  {
				   /*here puts integeral*/
				   F_AREA(A,f,tf);
				   /*  A_orintation=ND_DOT(A[0]/NV_MAG(A),A[1]/NV_MAG(A),A[2]/NV_MAG(A),0,0,1);/*to judge the orientation of A*/
				   d_vf=ND_DOT(A[0],A[1],A[2],C_U(c,t),C_V(c,t),C_W(c,t));/*dot product to calculate the volume flux*/
				    if (ND_DOT(C_U(c,t),C_V(c,t),C_W(c,t),0,0,1)<0)/*indicates area normal and vlocity orients the opposite, inflow*/
				     {
					   vf_in=vf_in+fabs(d_vf);
					   mf_in=mf_in+fabs(d_vf)*C_R(c,t);
					  /* fprintf(fp,"in %-5g\n",C_W(c,t) );  write the coordinate*/
					  q_in=q_in+C_CP(c,t)*C_R(c,t)*(C_T(c,t)-300)*fabs(d_vf);/*sum the heat flux*/					  
					  if ((fabs(xf[0]-0.3)<2e-1) && (epsilon <= fabs(C_P(c,t))))/*x= 0.1m-0.5m*/
					        {  
							   epsilon = fabs(C_P(c,t));
							   neutral_height=xf[1];
							}
					 }
					else 
						{ /*outflow*/
					    vf_out=vf_out+fabs(d_vf);
						mf_out=mf_out+fabs(d_vf)*C_R(c,t);
						/*fprintf(fp,"out %-5g\n",C_W(c,t) );  write the coordinate*/
						q_out=q_out+C_CP(c,t)*C_R(c,t)*(C_T(c,t)-300)*fabs(d_vf);/*sum the heat flux*/
						}
						cell_num++;
				  }
		       }
		   }
	    }
		end_c_loop(c,t)  
	}
	fprintf(fp,"%-5g %-5g %-5g %-5g %-5g %-5g %-5g %-5g %-5g\n",vf_out,mf_out,vf_in,mf_in,mf_in-mf_out,neutral_height,q_in,q_out,cell_num);  /*write the coordinate*/
	fclose(fp);
}



DEFINE_INIT(head_joint,d)
{
	Thread *t, *tf;
	cell_t c;
	face_t f;
	real xc[ND_ND];
	real xf[ND_ND];
	real A[ND_ND];  /*store face normal vector*/
	int n=-1,cell_num=0; /*loop variable*/
	
	FILE *fp; /*pointer to a file */
	fp = fopen("flux_joint.txt","w+");
	fprintf(fp,"/*this is recorded by FLUENT UDF*/\n");
	
	
	/*get the face normal area and record*/
	thread_loop_c(t,d)
	{
		begin_c_loop(c,t)
		{
		  C_CENTROID(xc,c,t);  /*get the center coordinates of one cell*/
		  if (xc[2]<=4.8 && xc[2]>=4.2 && xc[0]<=0.6 && xc[0]>=0.55) /*first, shrink search range*/
		  {
			   c_face_loop(c,t,n) /*loop over all the faces of one cell*/
			   {
			      f=C_FACE(c,t,n);
				  tf=C_FACE_THREAD(c,t,n);
				  F_CENTROID(xf,f,tf); /*get the coordinate of a face*/
				   
				 if ((fabs(xf[0]-0.6)<1e-3)&&(xc[0]<0.6)) /*indicates it is  boundary*/
				  {
				   /*here puts integeral*/
				   F_AREA(A,f,tf);
				   cell_num++;
				  }
		       }
		   }
	    }end_c_loop(c,t)
	}
	fprintf(fp,"face normal area is %g, %g, %g\t",A[0]/NV_MAG(A),A[1]/NV_MAG(A),A[2]/NV_MAG(A));
	 fprintf(fp,"cell_num is %d\t",cell_num);
	if (ND_DOT(A[0]/NV_MAG(A),A[1]/NV_MAG(A),A[2]/NV_MAG(A),1,0,0)>0)
	{
		fprintf(fp,"face area normal points to the branch\n");
	}
	else
	{
		fprintf(fp,"face area normal points to the main tunnel\n");
	}
	fprintf(fp,"vf_out mf_out vf_in mf_in mass_flow_in_rate neutral_height q_in q_out cell_num\n");
	fclose(fp);
}

DEFINE_EXECUTE_AT_END(flux_joint)
{
	Domain *d;
   	Thread *t,*tf;
    cell_t c;
	face_t f;
	real xc[ND_ND];
	real xf[ND_ND];
	real A[ND_ND];  /*store face normal vector*/

	real d_vf=0, d_mf=0;
	real vf_in=0, mf_in=0;
	real vf_out=0, mf_out=0;
	real q_in=0,q_out=0;
	real cell_num=0;
	
	int n=-1; /*loop variable*/
	real neutral_height=0,A_orintation=0,epsilon=100;
	
	FILE *fp; /*pointer to a file */
	fp = fopen("flux_joint.txt","a+");
	d = Get_Domain(1);
		
	thread_loop_c(t,d)
	{
		begin_c_loop(c,t)
		{
		  C_CENTROID(xc,c,t);  /*get the center coordinates of one cell*/
		 if (xc[2]<=4.8 && xc[2]>=4.2 && xc[0]<=0.6 && xc[0]>=0.55) /*first, shrink search range*/
		  {
			   c_face_loop(c,t,n) /*loop over all the faces of one cell*/
			   {
			      f=C_FACE(c,t,n);
				  tf=C_FACE_THREAD(c,t,n);
				  F_CENTROID(xf,f,tf); /*get the coordinate of a face*/
				  
				  if ((fabs(xf[0]-0.6)<1e-3)&&(xc[0]<0.6)) /*indicates it is  boundary*/
				  {
				   /*here puts integeral*/
				   F_AREA(A,f,tf);
				   /*  A_orintation=ND_DOT(A[0]/NV_MAG(A),A[1]/NV_MAG(A),A[2]/NV_MAG(A),0,0,1);/*to judge the orientation of A*/
				   d_vf=ND_DOT(A[0],A[1],A[2],C_U(c,t),C_V(c,t),C_W(c,t));/*dot product to calculate the volume flux*/
				    if (ND_DOT(C_U(c,t),C_V(c,t),C_W(c,t),1,0,0)<0)/*to main tunnel, inflow*/
				     {
					   vf_in=vf_in+fabs(d_vf);
					   mf_in=mf_in+fabs(d_vf)*C_R(c,t);
					  /* fprintf(fp,"in %-5g\n",C_W(c,t) );  write the coordinate*/
					  q_in=q_in+C_CP(c,t)*C_R(c,t)*(C_T(c,t)-300)*fabs(d_vf);/*sum the heat flux*/
					   if ((fabs(xf[2]-4.5)<2e-1) && (epsilon <= fabs(C_P(c,t))))/*x= 0.1m-0.5m*/
					        {  
							   epsilon = fabs(C_P(c,t));
							   neutral_height=xf[1];
							}
					 }
					else
						{ /*outflow*/
					    vf_out=vf_out+fabs(d_vf);
						mf_out=mf_out+fabs(d_vf)*C_R(c,t);
						/*fprintf(fp,"out %-5g\n",C_W(c,t) );  write the coordinate*/
						q_out=q_out+C_CP(c,t)*C_R(c,t)*(C_T(c,t)-300)*fabs(d_vf);/*sum the heat flux*/
						}
						cell_num++;
				  }
		       }
		   }
	    }end_c_loop(c,t)  
	}
	fprintf(fp,"%-5g %-5g %-5g %-5g %-5g %-5g %-5g %-5g %-5g\n",vf_out,mf_out,vf_in,mf_in,mf_in-mf_out,neutral_height,q_in,q_out,cell_num);  /*write the coordinate*/
	fclose(fp);
}

DEFINE_INIT(head_end,d)
{
	Thread *tf = Lookup_Thread(d,surfNO);/*get the branched end surface*/
	face_t f;
	real xf[ND_ND];
	real A[ND_ND];  /*store face normal vector*/
	int n=-1;        /*loop variable*/
	real face_num=0;
	
	FILE *fp; /*pointer to a file */
	fp = fopen("flux_end.txt","w+");
	fprintf(fp,"/*this is recorded by FLUENT UDF*/\n");
	
	begin_f_loop(f,tf)
	{
		F_CENTROID(xf,f,tf); /*get the coordinate of a face*/
		F_AREA(A,f,tf);      /*get the face normal vector*/
	    face_num++;
	}end_f_loop(f,tf)
	
	fprintf(fp,"face normal area is %g, %g, %g\t",A[0]/NV_MAG(A),A[1]/NV_MAG(A),A[2]/NV_MAG(A));
	fprintf(fp,"face_num is %g\t",face_num);
	
	if (ND_DOT(A[0]/NV_MAG(A),A[1]/NV_MAG(A),A[2]/NV_MAG(A),1,0,0)>0)
	{
		fprintf(fp,"face area normal points to the branch\n");
	}
	else
	{
		fprintf(fp,"face area normal points to the main tunnel\n");
	}
	fprintf(fp,"vf_out mf_out vf_in mf_in mass_flow_in_rate neutral_height q_in q_out face_num\n");
	fclose(fp);
}

DEFINE_EXECUTE_AT_END(flux_end)
{
	Domain *d=Get_Domain(1);
    Thread *tf = Lookup_Thread(d,surfNO);/*get the branched end surface*/
    Thread *t0, *t1=NULL;
	cell_t c0,c1=-1; /* adjacent cell */
	face_t f;
	real xf[ND_ND],xc[ND_ND];
	real A[ND_ND];  /*store face normal vector*/

	real d_vf=0, d_mf=0;
	real vf_in=0, mf_in=0;
	real vf_out=0, mf_out=0;
	real q_in=0,q_out=0;
	
	int n=-1; /*loop variable*/
	real neutral_height=0,A_orintation=0,epsilon=100,face_num=0;
	
	FILE *fp; /*pointer to a file */
	fp = fopen("flux_end.txt","a+");
	
	begin_f_loop(f,tf)
	{
		F_CENTROID(xf,f,tf); /*get the coordinate of a face*/
		F_AREA(A,f,tf);      /*get the face normal vector*/
	   
	   c0=F_C0(f,tf);
	   t0=F_C0_THREAD(f,tf);   
	   c1=F_C1(f,tf);
	   t1=F_C1_THREAD(f,tf);
	   /*fprintf(fp,"zone id= %-5g\n",Thread_ID(tf));  /*write the coordinate*/
	   C_CENTROID(xc,c0,t0);
	   
	   d_vf=ND_DOT(A[0],A[1],A[2],C_U(c0,t0),C_V(c0,t0),C_W(c0,t0));/*dot product to calculate the volume flux*/
	   if (ND_DOT(C_U(c0,t0),C_V(c0,t0),C_W(c0,t0),1,0,0)<0)/*to main tunnel, inflow*/
		   {
		      vf_in=vf_in+fabs(d_vf);
			  mf_in=mf_in+fabs(d_vf)*C_R(c0,t0);
			  /* fprintf(fp,"in %-5g\n",C_W(c,t) );  write the coordinate*/
			  q_in=q_in+C_CP(c0,t0)*C_R(c0,t0)*(C_T(c0,t0)-300)*fabs(d_vf);/*sum the heat flux*/
		    if (epsilon <= fabs(C_P(c0,t0)))
		       {
		  		   epsilon = fabs(C_P(c0,t0));
		  		   neutral_height=xf[1];
			   }
		   }
		else
			{ /*outflow*/
		      vf_out=vf_out+fabs(d_vf);
			  mf_out=mf_out+fabs(d_vf)*C_R(c0,t0);
			  /*fprintf(fp,"out %-5g\n",C_W(c,t) );  write the coordinate*/
			 q_out=q_out+C_CP(c0,t0)*C_R(c0,t0)*(C_T(c0,t0)-300)*fabs(d_vf);/*sum the heat flux*/
			}
			/*fprintf(fp,"x=%g,y=%g,z=%g\n",xc[0],xc[1],xc[2]);*/
			face_num++;
	}end_f_loop(f,tf)
	/*fprintf(fp,"ID=%d\n",THREAD_ID(tf));*/
	fprintf(fp,"%-5g %-5g %-5g %-5g %-5g %-5g %-5g %-5g %-5g\n",vf_out,mf_out,vf_in,mf_in,mf_in-mf_out,neutral_height,q_in,q_out,face_num);  /*write the coordinate*/
	fclose(fp);
}


/*---------------next record temperature gradient (smoke layer height) in main and branch---------------------------------------------*/

/*record layerheight, temperature gradient etc. */
DEFINE_INIT(layerheight_init,d)
{
	Thread *t, *tf;
	cell_t c;
	face_t f;
	real xc[ND_ND];
	real xf[ND_ND];
	real A[ND_ND];  /*store face normal vector*/
	int n=-1,cell_num=0; /*loop variable*/
	
	FILE *fp; /*pointer to a file */
	fp = fopen("smokelayerheight.txt","w+");
	fprintf(fp,"/*this is recorded by FLUENT UDF*/\n");
	fprintf(fp,"grad0,grad1,grad2,grad3,grad35,grad4,grad45,grad5,grad55,grad6,grad7,grad8,grad9,SLH0,SLH1,SLH2,SLH3,SLH35,SLH4,SLH45,SLH5,SLH55,SLH6,SLH7,SLH8,SLH9, vortex0,vortex1,vortex2,vortex3,vortex35,vortex4,vortex45,vortex5,vortex55,vortex6,vortex7,vortex8,vortex9");
	fprintf(fp,"\n");
	fclose(fp);
}

DEFINE_EXECUTE_AT_END(layerheight)
{
	Domain *d;
   	Thread *t,*tf;
    cell_t c;
	face_t f;
	real xc[ND_ND];
	real loopnum=0;
	real SLH0,SLH1,SLH2,SLH3,SLH35,SLH4,SLH45,SLH5,SLH55,SLH6,SLH7,SLH8,SLH9; /*SLH denotes smoke layer height*/
	real grad0,grad1,grad2,grad3,grad35,grad4,grad45,grad5,grad55,grad6,grad7,grad8,grad9;/*grad denote gradient*/
	real vortex0,vortex1,vortex2,vortex3,vortex35,vortex4,vortex45,vortex5,vortex55,vortex6,vortex7,vortex8,vortex9;/*vortex denote vortex*/
	FILE *fp; /*pointer to a file */
	
	SLH0=0,SLH1=0,SLH2=0,SLH3=0,SLH35=0,SLH4=0,SLH45=0,SLH5=0,SLH55=0,SLH6=0,SLH7=0,SLH8=0,SLH9=0;
	grad0=0,grad1=0,grad2=0,grad3=0,grad35=0,grad4=0,grad45=0,grad5=0,grad55=0,grad6=0,grad7=0,grad8=0,grad9=0;
	vortex0=0,vortex1=0,vortex2=0,vortex3=0,vortex35=0,vortex4=0,vortex45=0,vortex5=0,vortex55=0,vortex6=0,vortex7=0,vortex8=0,vortex9=0;
	

	fp = fopen("smokelayerheight.txt","a+");
	d = Get_Domain(1);
	
	thread_loop_c(t,d)
	{
		begin_c_loop(c,t)
		{
		  C_CENTROID(xc,c,t);  /*get the center coordinates of one cell*/
		  if (fabs(xc[0]-0.3)<=1e-2)  /*only consider the main tunnel central plane,the mesh resolution is 1e-2*/
		  {
			  if ((fabs(xc[2]-0)<=5e-2) && (grad0<C_T_G(c,t)[1]))  /*z=0 *//*gradient of temperature density, Y-component*/
			{
				 grad0=C_T_G(c,t)[1];
				 SLH0=xc[1];  /*get the height*/
				 vortex0=-C_DVDZ(c,t)+C_DWDY(c,t);/*calculate planner vorticity*/
			}
			else
			{
				if  ((fabs(xc[2]-1)<=5e-2) &&  (grad1<C_T_G(c,t)[1]))
				{
				  grad1=C_T_G(c,t)[1];
				  SLH1=xc[1];  /*get the height*/
				  vortex1=-C_DVDZ(c,t)+C_DWDY(c,t);/*calculate planner vorticity*/
				}
				else 
				{
					if ((fabs(xc[2]-2)<=5e-2) &&  (grad2<C_T_G(c,t)[1]))
					{
						grad2=C_T_G(c,t)[1];
				        SLH2=xc[1];  /*get the height*/
				        vortex2=-C_DVDZ(c,t)+C_DWDY(c,t);/*calculate planner vorticity*/
					}
					else
					{
						if ( (fabs(xc[2]-3)<=5e-2) && (grad3<C_T_G(c,t)[1]) )
						{
							grad3=C_T_G(c,t)[1];
				            SLH3=xc[1];  /*get the height*/
				            vortex3=-C_DVDZ(c,t)+C_DWDY(c,t);/*calculate planner vorticity*/
						}
						else
						{
							if ((fabs(xc[2]-3.5)<=5e-2)  && (grad35<C_T_G(c,t)[1]))
							{
								grad35=C_T_G(c,t)[1];
				                SLH35=xc[1];  /*get the height*/
				                vortex35=-C_DVDZ(c,t)+C_DWDY(c,t);/*calculate planner vorticity*/
							}
							else
							{
                               if ((fabs(xc[2]-4)<=5e-2) && (grad4<C_T_G(c,t)[1]))
							   {
								    grad4=C_T_G(c,t)[1];
				                    SLH4=xc[1];  /*get the height*/
				                    vortex4=-C_DVDZ(c,t)+C_DWDY(c,t);/*calculate planner vorticity*/   
							   }
							   else 
							   {
								   if ((fabs(xc[2]-4.5)<=5e-2) && (grad45<C_T_G(c,t)[1]))
								   {
									   grad45=C_T_G(c,t)[1];
				                       SLH45=xc[1];  /*get the height*/
				                       vortex45=-C_DVDZ(c,t)+C_DWDY(c,t);/*calculate planner vorticity*/ 
								   }
								   else
								   {
									   if ((fabs(xc[2]-5)<=5e-2) && (grad5<C_T_G(c,t)[1]))
									   {
										     grad5=C_T_G(c,t)[1];
				                             SLH5=xc[1];  /*get the height*/
				                             vortex5=-C_DVDZ(c,t)+C_DWDY(c,t);/*calculate planner vorticity*/										   
									   }
									   else
									   {
										   if ((fabs(xc[2]-5.5)<=5e-2) && (grad55<C_T_G(c,t)[1]))
										   {
											   grad55=C_T_G(c,t)[1];
				                               SLH55=xc[1];  /*get the height*/
				                                vortex55=-C_DVDZ(c,t)+C_DWDY(c,t);/*calculate planner vorticity*/
										   }
										   else
										   {
											   if ((fabs(xc[2]-6)<=5e-2) && (grad6<C_T_G(c,t)[1]))
											   {
												    grad6=C_T_G(c,t)[1];
				                                    SLH6=xc[1];  /*get the height*/
				                                    vortex6=-C_DVDZ(c,t)+C_DWDY(c,t);/*calculate planner vorticity*/
											   }
											   else
											   {
												   if ((fabs(xc[2]-7)<=5e-2) && (grad7<C_T_G(c,t)[1]))
												   {
													   grad7=C_T_G(c,t)[1];
				                                       SLH7=xc[1];  /*get the height*/
				                                       vortex7=-C_DVDZ(c,t)+C_DWDY(c,t);/*calculate planner vorticity*/
												   }
												   else
												   {
													   if ((fabs(xc[2]-8)<=5e-2)  &&  (grad8<C_T_G(c,t)[1]))
													   {
														   grad8=C_T_G(c,t)[1];
				                                           SLH8=xc[1];  /*get the height*/
				                                            vortex8=-C_DVDZ(c,t)+C_DWDY(c,t);/*calculate planner vorticity*/
													   }
													   else 
													   {
														   if ((fabs(xc[2]-9)<=5e-2) && (grad9<C_T_G(c,t)[1]) )
														   {
															    grad9=C_T_G(c,t)[1];
				                                                SLH9=xc[1];  /*get the height*/
			                                                    vortex9=-C_DVDZ(c,t)+C_DWDY(c,t);/*calculate planner vorticity*/
														   }
													   }
													   
												   }
											   }
										   }
										   
									   }
								   }
							   }
							   
							}
							
						}
					}
				}
			}
		  
		 }
	    }
		end_c_loop(c,t)  
	}
	fprintf(fp,"%-5g %-5g %-5g %-5g %-5g %-5g %-5g %-5g %-5g %-5g %-5g %-5g %-5g %-5g %-5g %-5g %-5g %-5g %-5g %-5g %-5g %-5g %-5g %-5g %-5g %-5g %-5g %-5g %-5g %-5g %-5g %-5g %-5g %-5g %-5g %-5g %-5g %-5g %-5g\n", grad0,grad1,grad2,grad3,grad35,grad4,grad45,grad5,grad55,grad6,grad7,grad8,grad9,SLH0,SLH1,SLH2,SLH3,SLH35,SLH4,SLH45,SLH5,SLH55,SLH6,SLH7,SLH8,SLH9, vortex0,vortex1,vortex2,vortex3,vortex35,vortex4,vortex45,vortex5,vortex55,vortex6,vortex7,vortex8,vortex9);
	fclose(fp);
}

DEFINE_INIT(layerheight_branch_init,d)
{
	Thread *t, *tf;
	cell_t c;
	face_t f;
	real xc[ND_ND];
	real xf[ND_ND];
	real A[ND_ND];  /*store face normal vector*/
	int n=-1,cell_num=0; /*loop variable*/
	
	FILE *fp; /*pointer to a file */
	fp = fopen("smokelayerheight_branch.txt","w+");
	fprintf(fp,"/*this is recorded by FLUENT UDF*/\n");
	fprintf(fp,"grad06,grad16,grad26,grad36,grad46,SLH06,SLH16,SLH26,SLH36,SLH46,vortex06,vortex16,vortex26,vortex36,vortex46");
	fprintf(fp,"\n");
	fclose(fp);
}

DEFINE_EXECUTE_AT_END(layerheight_branch)
{
	Domain *d;
   	Thread *t,*tf;
    cell_t c;
	face_t f;
	real xc[ND_ND];
	real SLH06,SLH16,SLH26,SLH36,SLH46; /*SLH denotes smoke layer height*/
	real grad06,grad16,grad26,grad36,grad46;/*grad denote gradient*/
	real vortex06,vortex16,vortex26,vortex36,vortex46;/*vortex denote vortex*/
	FILE *fp; /*pointer to a file */
	
	SLH06=0,SLH16=0,SLH26=0,SLH36=0,SLH46=0;
	grad06=0,grad16=0,grad26=0,grad36=0,grad46=0;
	vortex06=0,vortex16=0,vortex26=0,vortex36=0,vortex46=0;
	

	fp = fopen("smokelayerheight_branch.txt","a+");
	d = Get_Domain(1);
	
	thread_loop_c(t,d)
	{
		begin_c_loop(c,t)
		{
		  C_CENTROID(xc,c,t);  /*get the center coordinates of one cell*/
		  if (fabs(xc[2]-4.5)<=1e-2)  /*only consider the central plane*/
		  {
			if ((fabs(xc[0]-0.6)<=5e-2)  &&  (grad06<C_T_G(c,t)[1]))  /*x=0.6 */
			{
				  grad06=C_T_G(c,t)[1];
				  SLH06=xc[1];  /*get the height*/
				  vortex06=-C_DVDZ(c,t)+C_DWDY(c,t);/*calculate planner vorticity*/
			}
			else
			{
				if ((fabs(xc[0]-1.6)<=5e-2) && (grad16<C_T_G(c,t)[1]))  /*x=1.6 */
			    {
				  grad16=C_T_G(c,t)[1];
				  SLH16=xc[1];  /*get the height*/
				  vortex16=-C_DVDZ(c,t)+C_DWDY(c,t);/*calculate planner vorticity*/
			    }
				else
				{
					if ((fabs(xc[0]-2.6)<=5e-2)  && (grad26<C_T_G(c,t)[1]))  /*x=2.6 */
			       {
				     grad26=C_T_G(c,t)[1];
				     SLH26=xc[1];  /*get the height*/
				     vortex26=-C_DVDZ(c,t)+C_DWDY(c,t);/*calculate planner vorticity*/
				   }
				   else 
				   {
					   if ((fabs(xc[0]-3.6)<=5e-2) && (grad36<C_T_G(c,t)[1])) /*x=3.6 */
					   {
						   grad36=C_T_G(c,t)[1];
				           SLH36=xc[1];  /*get the height*/
				           vortex36=-C_DVDZ(c,t)+C_DWDY(c,t);/*calculate planner vorticity*/ 
					   }
					   else
					   {
						  if ((fabs(xc[0]-4.6)<=5e-2) && (grad46<C_T_G(c,t)[1])) /*x=4.6 */
						  {
							  grad46=C_T_G(c,t)[1];
				              SLH46=xc[1];  /*get the height*/
				              vortex46=-C_DVDZ(c,t)+C_DWDY(c,t);/*calculate planner vorticity*/
						  }
					   }
				   }
			    }
			}
		 }
		}end_c_loop(c,t)  	
	}
	
  fprintf(fp,"%-5g %-5g %-5g %-5g %-5g %-5g %-5g %-5g %-5g %-5g %-5g %-5g %-5g %-5g %-5g\n", grad06,grad16,grad26,grad36,grad46,SLH06,SLH16,SLH26,SLH36,SLH46,vortex06,vortex16,vortex26,vortex36,vortex46);
  fclose(fp);	
}

/*report the inertial force */
/*first write txt head*/
DEFINE_INIT(txt_inertial_force,d)
{
	Thread *t, *tf;
	cell_t c;
	face_t f;
	real xc[ND_ND];
	real xf[ND_ND];
	real A[ND_ND];  /*store face normal vector*/
	int n=-1,cell_num=0; /*loop variable*/
	
	FILE *fp; /*pointer to a file */
	fp = fopen("inertial_force.txt","w+");
	
	fprintf(fp,"P_xmax P_xmin F_Y P_zmax P_zmin cell_num_xmax  cell_num_xmin  cell_num_zmax  cell_num_zmin\n");
	fclose(fp);
}

DEFINE_EXECUTE_AT_END(record_inertial_force)
 {
    Domain *d;
   	Thread *t,*tf;
    cell_t c;
	face_t f;
	real xc[ND_ND],xf[ND_ND];
	real A[ND_ND],cell_num_xmax=0,cell_num_xmin=0,cell_num_zmax=0,cell_num_zmin=0;

	real D,lf,xmax,xmin,zmin,zmax,source,volume,Q_xing;
	real DP_xmax=0,DP_xmin=0,DP_zmax=0,DP_zmin=0;
	real P_xmax=0,P_xmin=0,P_zmax=0,P_zmin=0,F_Y=0,f_buoyancy=0;
	real height_crit=0;
	int n=-1;
	
	FILE *fp; /*pointer to a file */
	fp = fopen("inertial_force.txt","a+");
	d = Get_Domain(1);       /* Get the domain using ANSYS FLUENT utility */
	
    D=2*(a+b)/M_PI;   /*D is perimeter diamter*/
	Q_xing=Q/(1.2*1.0*300*sqrt(9.81)*pow(D,2.5));
	lf=(3.7*pow(Q_xing,0.4)-1.02)*D;  /*calculate the height of fire*/
	if(lf>0.45)
	{	lf=0.45;}
	
	volume=a*b*lf;/*calculate the volume of fire source*/
	
    xmax=0.3+a/2;
    xmin=0.3-a/2;
    zmax=4.5+b/2;
    zmin=4.5-b/2;
	height_crit=0.5*lf;

   /* Loop over all cell threads in the domain */
   thread_loop_c(t,d)
      {
      /* Loop over all cells */
      begin_c_loop(c,t)
        {
         C_CENTROID(xc,c,t);  /*get the center coordinates of one cell*/
		 if (xc[0]<=xmax && xc[0]>=xmin && xc[2]<=zmax && xc[2]>=zmin && xc[1]<=height_crit)  /*whether flame region*/
		  {
			   f_buoyancy=(1.1774-C_R(c,t))*9.81*C_VOLUME(c,t); /*calculate the buoyancy force of fire grids*/
			   
			    c_face_loop(c,t,n) /*loop over all the faces of one cell*/
			   {
			      f=C_FACE(c,t,n);
				  tf=C_FACE_THREAD(c,t,n);
				  F_CENTROID(xf,f,tf); /*get the coordinate of a face*/
				  
				 if ((fabs(xf[0]-(xmax))<1e-4))/*indicates it is  boundary*/
				  {
					 cell_num_xmax++;
				   /*here puts integeral*/
				   F_AREA(A,f,tf);
				   DP_xmax=NV_MAG(A)*(F_P(f,tf)+0.5*C_R(c,t)*SQR(C_U(c,t)));/*pressure of a small element face*/	
				  }  
				else
					{
						if ((fabs(xf[0]-(xmin))<1e-4))/*indicates it is  boundary*/
						{
						  cell_num_xmin++;
					     /*here puts integeral*/
				         F_AREA(A,f,tf);
				          DP_xmin=NV_MAG(A)*(F_P(f,tf)+0.5*C_R(c,t)*SQR(C_U(c,t)));/*pressure of a small element face*/		
					    }
					   else
					   {
						  if  ((fabs(xf[2]-(zmax))<1e-4))/*indicates it is  boundary*/
						  {
							  cell_num_zmax++;
							  F_AREA(A,f,tf);
				             DP_zmax=NV_MAG(A)*(F_P(f,tf)+0.5*C_R(c,t)*SQR(C_W(c,t)));/*pressure of a small element face*/		
						  }
						  else
						  {
							  if ((fabs(xf[2]-(zmin))<1e-4))/*indicates it is  boundary*/
							  {
								  cell_num_zmin++;
								  F_AREA(A,f,tf);
				                  DP_zmin=NV_MAG(A)*(F_P(f,tf)+0.5*C_R(c,t)*SQR(C_W(c,t)));/*pressure of a small element face*/		
							  }
						  }
					  }
				   }
			   }
		 }   
		  P_xmax=P_xmax+DP_xmax;
		  P_xmin=P_xmin+DP_xmin;
	      F_Y=F_Y+f_buoyancy;
		  P_zmax=P_zmax+DP_zmax;
		  P_zmin=P_zmin+DP_zmin;
        }end_c_loop(c,t)
    fprintf(fp,"%-5g %-5g %-5g %-5g %-5g %-5g %-5g %-5g %-5g\n",P_xmax,P_xmin,F_Y,P_zmax,P_zmin,cell_num_xmax,cell_num_xmin, cell_num_zmax, cell_num_zmin);  /*write the coordinate*/
	fclose(fp);
	}
}


DEFINE_ON_DEMAND(report_inertial_force)
{
	/*this part may be not right*/
	Domain *d;
   	Thread *t;
    cell_t c;
	real xc[ND_ND],cell_num=0;

	real f_x=0,f_y=0,f_z=0;
	real f_x_sum=0,f_y_sum=0,f_z_sum=0,f_sum=0;
	real D,lf,xmax,xmin,zmin,zmax,source,volume,Q_xing;
	real DrhoDX=0,DrhoDY=0,DrhoDZ=0;
	real transient_term_x=0,transient_term_y=0,transient_term_z=0;
	
	FILE *fp; /*pointer to a file */
	fp = fopen("inertial_force.txt","a+");
	d = Get_Domain(1);
	
    D=2*(a+b)/M_PI;   /*D is perimeter diamter*/
	Q_xing=Q/(1.2*1.0*300*sqrt(9.81)*pow(D,2.5));
	lf=(3.7*pow(Q_xing,0.4)-1.02)*D;  /*calculate the height of fire*/
	if(lf>0.45)
	{	lf=0.45;}

	volume=a*b*lf;/*calculate the volume of fire source*/
	
    xmax=0.3+a/2;
    xmin=0.3-a/2;
    zmax=4.5+b/2;
    zmin=4.5-b/2;
    
	thread_loop_c(t,d)
	{
		begin_c_loop(c,t)
		{
		  C_CENTROID(xc,c,t);  /*get the center coordinates of one cell*/
		 if (xc[0]<=xmax && xc[0]>=xmin && xc[2]<=zmax && xc[2]>=zmin && xc[1]<=lf)  /*whether flame region*/
		  {
			   DrhoDX=(C_P_G(c,t)[0]-C_P(c,t)/C_T(c,t)*C_T_G(c,t)[0])/Rg/C_T(c,t);
			   DrhoDY=(C_P_G(c,t)[1]-C_P(c,t)/C_T(c,t)*C_T_G(c,t)[1])/Rg/C_T(c,t);
			   DrhoDZ=(C_P_G(c,t)[2]-C_P(c,t)/C_T(c,t)*C_T_G(c,t)[2])/Rg/C_T(c,t);
			   cell_num++;
			   transient_term_x=(C_R(c,t)*C_U(c,t)-C_U_M1(c,t)*C_R_M1(c,t))/CURRENT_TIMESTEP;
			   transient_term_y=(C_R(c,t)*C_V(c,t)-C_V_M1(c,t)*C_R_M1(c,t))/CURRENT_TIMESTEP;
			   transient_term_z=(C_R(c,t)*C_W(c,t)-C_W_M1(c,t)-C_R_M1(c,t))/CURRENT_TIMESTEP;
			   f_x=C_VOLUME(c,t)*(transient_term_x+C_U(c,t)*C_U(c,t)*DrhoDX+C_U(c,t)*C_V(c,t)*DrhoDY+C_U(c,t)*C_W(c,t)*DrhoDZ+2*C_R(c,t)*C_U(c,t)*C_DUDX(c,t)+C_R(c,t)*C_U(c,t)*C_DVDY(c,t)+C_R(c,t)*C_U(c,t)*C_DWDZ(c,t)+C_R(c,t)*C_V(c,t)*C_DUDY(c,t)+C_R(c,t)*C_W(c,t)*C_DUDZ(c,t));
			   f_y=C_VOLUME(c,t)*(transient_term_y+C_V(c,t)*C_U(c,t)*DrhoDX+C_V(c,t)*C_V(c,t)*DrhoDY+C_V(c,t)*C_W(c,t)*DrhoDZ+2*C_R(c,t)*C_V(c,t)*C_DVDY(c,t)+C_R(c,t)*C_U(c,t)*C_DVDX(c,t)+C_R(c,t)*C_V(c,t)*C_DUDX(c,t)+C_R(c,t)*C_W(c,t)*C_DVDZ(c,t)+C_R(c,t)*C_V(c,t)*C_DWDZ(c,t));
			   f_z=C_VOLUME(c,t)*(transient_term_z+C_W(c,t)*C_U(c,t)*DrhoDX+C_W(c,t)*C_V(c,t)*DrhoDY+C_W(c,t)*C_W(c,t)*DrhoDZ+2*C_R(c,t)*C_W(c,t)*C_DWDZ(c,t)+C_R(c,t)*C_W(c,t)*C_DUDX(c,t)+C_R(c,t)*C_U(c,t)*C_DWDX(c,t)+C_R(c,t)*C_W(c,t)*C_DVDY(c,t)+C_R(c,t)*C_V(c,t)*C_DWDY(c,t));
	      }
		  f_x_sum=f_x_sum+f_x;
		  f_y_sum=f_y_sum+f_y;
		  f_z_sum=f_z_sum+f_z;
		  f_sum=sqrt(SQR(f_x_sum)+SQR(f_y_sum)+SQR(f_z_sum));
		}end_c_loop(c,t)
	}
	fprintf(fp,"%-5g %-5g %-5g %-5g %-5g %-5g\n",f_x_sum,f_y_sum,f_z_sum,f_sum,volume,lf);  /*write the coordinate*/
	fclose(fp);
}
DEFINE_INIT(init,d)
{
    real xc[ND_ND];
    real D,Q_xing,lf,volume,heat_source,mass_source,source2;
 
    D=2*(a+b)/M_PI;   /*D is perimeter diamter*/
	Q_xing=Q/(1.2*1.0*300*sqrt(9.81)*pow(D,2.5));
	lf=(3.7*pow(Q_xing,0.4)-1.02)*D;  /*calculate the height of fire*/
	lf=(lf>0.45)?0.45:lf;
	volume=a*b*lf;/*calculate the volume of fire source*/
	
    heat_source=Q*1000/volume;
	
		if ((3.7*pow(Q_xing,0.4)-1.02)*D>0.45) /*flame touch the ceiling*/
		    {mass_source=0.284*pow(Q,1.0/3)*pow(lf,2.0/3)/a/pow(b,3.0)*pow(xc[2]-4.5,2.0)+0.0473*pow(Q,1.0/3)*pow(lf,2.0/3)/(a*b);}
		else /*flame not touch the ceiling*/
		{mass_source=0.071*pow(Q,1.0/3)*pow(lf,5.0/3)/volume;}
		  
	CX_Message("mass_source=%-5g\n",mass_source);
	CX_Message("(3.7*pow(Q_xing,0.4)-1.02)=%-5g\n",(3.7*pow(Q_xing,0.4)-1.02));
	CX_Message("heat_source=%-5g\n",heat_source);
}