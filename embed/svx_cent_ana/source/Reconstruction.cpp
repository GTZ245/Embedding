#include "Reconstruction.h"
#include "MyEvent.h"
#include <TFile.h>
#include <TNtuple.h>
#include <TMatrixD.h>
#include <iostream>
#include <TMath.h>

using namespace std;

namespace DileptonAnalysis
{

struct PointValue_3x3
{
	float x,y,z;// independent variables
	float dep[3];// dependent variables
	PointValue_3x3(float X, float Y, float Z, float A, float B, float C) : x(X), y(Y), z(Z)
	{
		dep[0] = A;
		dep[1] = B;
		dep[2] = C;
	}
	PointValue_3x3(){}
	PointValue_3x3( PointValue_3x3 const& other )
	{
		x = other.x;
		y = other.y;
		z = other.z;
		dep[0] = other.dep[0];
		dep[1] = other.dep[1];
		dep[2] = other.dep[2];
	}
	PointValue_3x3& operator= (const PointValue_3x3& other)
        {
		x = other.x;
		y = other.y;
		z = other.z;
		dep[0] = other.dep[0];
		dep[1] = other.dep[1];
		dep[2] = other.dep[2];
		return *this;
         }
};


class PointVal3x3Sorted
{
  protected:
    float x_lo, x_hi;
    float y_lo, y_hi;
    float z_lo, z_hi;

    int level;
    int maxlevel;

    PointVal3x3Sorted* containers[2][2][2];

    public:
      PointVal3x3Sorted( float X_LO, float X_HI, float Y_LO, float Y_HI, float Z_LO, float Z_HI, int MLEV, int LEV=0 ) : x_lo(X_LO), x_hi(X_HI), y_lo(Y_LO), y_hi(Y_HI), z_lo(Z_LO), z_hi(Z_HI), level(LEV), maxlevel(MLEV)
      {
        for(unsigned int i=0;i<2;++i){for(unsigned int j=0;j<2;++j){for(unsigned int k=0;k<2;++k){containers[i][j][k]=nullptr;}}}
      }
      virtual ~PointVal3x3Sorted()
      {
        for(unsigned int i=0;i<2;++i){
          for(unsigned int j=0;j<2;++j){
          	for(unsigned int k=0;k<2;++k){
              if(containers[i][j][k]!=nullptr)
              {
                delete containers[i][j][k];
              }}}}
      }
      virtual bool insert( PointValue_3x3 const& point );
      virtual void append_list( vector<PointValue_3x3>& point_list, float X_LO, float X_HI, float Y_LO, float Y_HI, float Z_LO, float Z_HI )
      {
        for(unsigned int i=0;i<2;++i){
          for(unsigned int j=0;j<2;++j){
          	for(unsigned int k=0;k<2;++k){
              if(containers[i][j][k]==nullptr){continue;}
              if( (containers[i][j][k]->x_hi<X_LO) || (containers[i][j][k]->x_lo>X_HI) || (containers[i][j][k]->y_hi<Y_LO) || (containers[i][j][k]->y_lo>Y_HI) || (containers[i][j][k]->z_hi<Z_LO) || (containers[i][j][k]->z_lo>Z_HI) ){continue;}
              containers[i][j][k]->append_list( point_list, X_LO, X_HI, Y_LO, Y_HI, Z_LO, Z_HI );}}}
      }
};

class PointVal3x3SortedEnd : public PointVal3x3Sorted
{
  public:
    PointVal3x3SortedEnd( float X_LO, float X_HI, float Y_LO, float Y_HI, float Z_LO, float Z_HI, int MLEV, int LEV=0 ) : PointVal3x3Sorted( X_LO,X_HI,Y_LO,Y_HI,Z_LO,Z_HI,MLEV,LEV ){}
    virtual ~PointVal3x3SortedEnd(){ points.clear(); }
    virtual bool insert( PointValue_3x3 const& point )
    {
      points.push_back(point);
      return true;
    }
    virtual void append_list( vector<PointValue_3x3>& point_list, float X_LO, float X_HI, float Y_LO, float Y_HI, float Z_LO, float Z_HI )
    {
      for(unsigned int i=0;i<points.size();++i)
      {
        if( (points[i].x<X_LO) || (points[i].x>X_HI) || (points[i].y<Y_LO) || (points[i].y>Y_HI) || (points[i].z<Z_LO) || (points[i].z>Z_HI) ){continue;}
        point_list.push_back( points[i] );
      }
    }
  protected:
    vector<PointValue_3x3> points;
};

bool PointVal3x3Sorted::insert( PointValue_3x3 const& point )
{
  if( (point.x < x_lo) || (point.y < y_lo) || (point.z < z_lo) || (point.x > x_hi) || (point.y > y_hi) || (point.z > z_hi) )
  {
    return false;
  }

  int x_ind = 0;
  if(point.x > (x_lo + 0.5*(x_hi-x_lo))){x_ind=1;}
  int y_ind = 0;
  if(point.y > (y_lo + 0.5*(y_hi-y_lo))){y_ind=1;}
  int z_ind = 0;
  if(point.z > (z_lo + 0.5*(z_hi-z_lo))){z_ind=1;}

  if( containers[x_ind][y_ind][z_ind] == nullptr )
  {
    float x_lo_new = x_lo + (float(x_ind))*0.5*(x_hi-x_lo);
    float x_hi_new = x_lo_new + 0.5*(x_hi-x_lo);

    float y_lo_new = y_lo + (float(y_ind))*0.5*(y_hi-y_lo);
    float y_hi_new = y_lo_new + 0.5*(y_hi-y_lo);

    float z_lo_new = z_lo + (float(z_ind))*0.5*(z_hi-z_lo);
    float z_hi_new = z_lo_new + 0.5*(z_hi-z_lo);

    if(level < maxlevel)
    {
      containers[x_ind][y_ind][z_ind] = new PointVal3x3Sorted( x_lo_new,x_hi_new, y_lo_new,y_hi_new, z_lo_new,z_hi_new, maxlevel, level+1 );
    }
    else
    {
      containers[x_ind][y_ind][z_ind] = new PointVal3x3SortedEnd( x_lo_new,x_hi_new, y_lo_new,y_hi_new, z_lo_new,z_hi_new, maxlevel, level+1 );
    }
  }
  return containers[x_ind][y_ind][z_ind]->insert( point );
}


class Reconstruction::impl
{
	public:
		impl( const char* lookup_file_name ) : alpha_r_z( 0., 0.7, 0., 150., 0., 100., 6 )
		{
			TFile f(lookup_file_name);
			TNtuple* t=0;
			f.GetObject("alpha_p_r_phi_theta_z", t);
			float alpha, p, r, phi, theta, z;
			t->SetBranchAddress("alpha",&alpha);
			t->SetBranchAddress("p",&p);
			t->SetBranchAddress("r",&r);
			t->SetBranchAddress("phi",&phi);
			t->SetBranchAddress("theta",&theta);
			t->SetBranchAddress("z",&z);
			for(int i=0;i<t->GetEntries();i+=1)
			{
				t->GetEntry(i);

				if ( z < 0.) continue;
				if ( alpha > 0. ) continue;
				if ( phi < 0. || phi > TMath::Pi() ) continue;
				if ( theta < 0. ) continue;

				PointValue_3x3 point(TMath::Abs(alpha),r,TMath::Abs(z), TMath::Abs(phi), p, TMath::Abs(theta));
				alpha_r_z.insert(point);
			}
			f.Close();
		}

		PointVal3x3Sorted alpha_r_z;
};

Reconstruction::Reconstruction( const char* lookup_file_name )
{
	pimpl = new impl(lookup_file_name);
}

Reconstruction::~Reconstruction()
{
	delete pimpl;
}




static const float max_r = 30.;
static const float lookup_alpha_delta_base = 0.01;
static const float lookup_r_delta_base = 2.0;
static const float lookup_z_delta_base = 4.0;

static float evaluate_fit( TMatrixD& beta, float alpha, float r, float z_DC_m_z_ver )
{
	float retval = 0.;
	retval += beta(0,0);
	retval += beta(1,0)*alpha;
	retval += beta(2,0)*r;
	retval += beta(3,0)*z_DC_m_z_ver;
	retval += beta(4,0)*alpha*alpha;
	retval += beta(5,0)*alpha*r;
	retval += beta(6,0)*alpha*z_DC_m_z_ver;
	retval += beta(7,0)*r*r;
	retval += beta(8,0)*r*z_DC_m_z_ver;
	retval += beta(9,0)*z_DC_m_z_ver*z_DC_m_z_ver;

	return retval;
}

static void remove_outliers( int ind, TMatrixD& beta, vector<PointValue_3x3>& points_temp, vector<PointValue_3x3>& points )
{
	points.clear();

	float sig = 0.;
	for( unsigned int i=0; i<points_temp.size(); ++i )
	{
		float fval = evaluate_fit( beta, points_temp[i].x, points_temp[i].y, points_temp[i].z );
		sig += ( fval - points_temp[i].dep[ind] )*( fval - points_temp[i].dep[ind] );
	}
	sig /= ((float)(points_temp.size()));
	sig = sqrt(sig);

	for( unsigned int i=0; i<points_temp.size(); ++i )
	{
		float fval = evaluate_fit( beta, points_temp[i].x, points_temp[i].y, points_temp[i].z );
		if( fabs( fval - points_temp[i].dep[ind] ) < 3.*sig ){points.push_back(points_temp[i]);}
	}
}

TMatrixD perform_fit( int ind, vector<PointValue_3x3> const& points )
{
	TMatrixD X( points.size(), 10 );
	TMatrixD y( points.size(), 1  );
	for(unsigned int i=0;i<points.size();i+=1)
	{
		y(i,0) = points[i].dep[ind];
		X(i,0) = 1.;
		X(i,1) = points[i].x;
		X(i,2) = points[i].y;
		X(i,3) = points[i].z;
		X(i,4) = points[i].x * points[i].x;
		X(i,5) = points[i].x * points[i].y;
		X(i,6) = points[i].x * points[i].z;
		X(i,7) = points[i].y * points[i].y;
		X(i,8) = points[i].y * points[i].z;
		X(i,9) = points[i].z * points[i].z;
	}

	TMatrixD Xt( 10, points.size() );
	Xt.Transpose(X);

	TMatrixD XtX = Xt * X;
	XtX.Invert();

	TMatrixD beta = XtX * (Xt * y);

	return beta;
}


// 0 <--> phi
// 1 <--> p
// 2 <--> theta
static float lookup_fit( float alpha, float r, float z_DC_m_z_ver, int ind, PointVal3x3Sorted& alpha_r_z, float& lookup_ind )
{
	// r can be -9999 when find_conv() fails
	if ( r<0 || r>max_r ) { lookup_ind = -9999.; return 0; }
	if ( alpha<0 || z_DC_m_z_ver<0 )
	{
		cout<<"input error!!!"<<endl;
		cout<<"lookup table contain only e- entries, please convert e+ to e- before fitting with lookup table"<<endl;
		lookup_ind = -9999;
		return 0;
	}

	float lookup_alpha_delta = lookup_alpha_delta_base;
	float lookup_r_delta = lookup_r_delta_base;
	float lookup_z_delta = lookup_z_delta_base;

	vector<PointValue_3x3> points;

	alpha_r_z.append_list( points, alpha - lookup_alpha_delta, alpha + lookup_alpha_delta, r - lookup_r_delta, r + lookup_r_delta, z_DC_m_z_ver - lookup_z_delta, z_DC_m_z_ver + lookup_z_delta );

	if( points.size() < 20 )
	{
		lookup_ind = -9999.;
		return false;
	}

	// p0 + p1*x + p2*y + p3*z + p4*x^2 + p5*xy + p6*xz + p7*y^2 + p8*yz + p9*z^2

	TMatrixD beta = perform_fit( ind, points );
	{
		vector<PointValue_3x3> points_temp = points;
		remove_outliers( ind, beta, points_temp, points );
	}

	if( points.size() < 20 )
	{
		lookup_ind = -9999.;
		return false;
	}
	
	TMatrixD beta2 = perform_fit( ind, points );

	lookup_ind = evaluate_fit( beta2, alpha, r, z_DC_m_z_ver );
	return true;
}


static bool project_phi( float alpha, float phiDC, float r, float z_DC_m_z_ver, PointVal3x3Sorted& alpha_r_z, float& phi_r )
{
	//DC coord sys: phi -0.5Pi~1.5Pi
	if ( (r<0) || (r>max_r) ) { phi_r=-9999.; return 0; }

	float lookup_phi = 0;
	if ( !lookup_fit(TMath::Abs(alpha), r, TMath::Abs(z_DC_m_z_ver), 0, alpha_r_z, lookup_phi) ) { phi_r=-9999.; return 0; }
	
	if ( alpha>0 )
	{ // electron!
		if ( lookup_phi<0 )
		{
			cout<<"error!!! phi_conv-phi_DC can't be negative for electron using ++ field"<<endl;
			cout<<"lookup_phi "<<lookup_phi<<endl;
			cout<<"alpha "<<alpha<<endl;
			cout<<"r "<<r<<endl;
			cout<<"z_DC_m_z_ver "<<TMath::Abs(z_DC_m_z_ver)<<endl;
			phi_r = -9999.;
			return 0;
		} 	
	}
	else 
	{ // positron!
		lookup_phi = - lookup_phi;
		if ( lookup_phi>0 )
		{
			cout<<"error!!! phi_conv-phi_DC can't be positive for positron using ++ field"<<endl;
			cout<<"lookup_phi "<<lookup_phi<<endl;
			cout<<"alpha "<<alpha<<endl;
			cout<<"r "<<r<<endl;
			cout<<"z_DC_m_z_ver "<<TMath::Abs(z_DC_m_z_ver)<<endl;
			phi_r = -9999.;
			return 0;
		} 
	}

	phi_r = phiDC + lookup_phi;
	// periodic boundary, to ensure phi_r from [-0.5pi, 1.5pi]
	if (phi_r>1.5*TMath::Pi())  phi_r = phi_r-2*TMath::Pi();
	if (phi_r<-0.5*TMath::Pi()) phi_r = phi_r+2*TMath::Pi();
	// else already in range, phi_r = phi_r
	return 1;
}

static bool project_theta( float alpha, float r, float z_DC_m_z_ver, PointVal3x3Sorted& alpha_r_z, float& theta_r )
{
	//DC coord sys: phi -0.5Pi~1.5Pi
	if ( (r<0) || (r>max_r) ) { theta_r=-9999.; return 0; }

	float lookup_theta = 0;
	if ( !lookup_fit(TMath::Abs(alpha), r, TMath::Abs(z_DC_m_z_ver), 2, alpha_r_z, lookup_theta) ) { theta_r = -9999.; return 0; }
	// cout<<"TMath::Abs(z_DC_m_z_ver) "<<TMath::Abs(z_DC_m_z_ver)<<" lookup_theta "<<lookup_theta<<endl;

	if ( lookup_theta<0 )
	{
		cout<<"error!!! there should be no negative theta_conv from lookup table"<<endl;
		theta_r = -9999.;
		return 0;
	} 
	if ( z_DC_m_z_ver>0 ) theta_r = lookup_theta; // top half
	if ( z_DC_m_z_ver<0 ) theta_r = TMath::Pi() - lookup_theta; // bottom half
	return 1;
}

static bool delta_phi( float alpha_e, float alpha_p, float phi_e, float phi_p, float r, float z_DC_m_z_ver_e, float z_DC_m_z_ver_p, PointVal3x3Sorted& alpha_r_z, float& dphi )
{
	float phir_e, phir_p;
	if ( !project_phi(alpha_e, phi_e, r, z_DC_m_z_ver_e, alpha_r_z, phir_e) || !project_phi(alpha_p, phi_p, r, z_DC_m_z_ver_p, alpha_r_z, phir_p) ) { dphi = -9999.; return 0; }

	dphi = phir_p - phir_e;
	// the dphi we use here should alway be |dphi|<=pi
	if ( dphi>TMath::Pi() )  dphi = 2*TMath::Pi()-dphi; // 1.3pi->0.7pi
	if ( dphi<-TMath::Pi() ) dphi = -2*TMath::Pi()-dphi; // -1.3pi->-0.7pi
	return 1;
}

static bool find_rconv( float alpha_e, float alpha_p, float phi_e, float phi_p, float z_DC_m_z_ver_e, float z_DC_m_z_ver_p, PointVal3x3Sorted& alpha_r_z, float& r_conv )
{ // assumption: delta_phi(r) should be monotomic increase within [0, max_r]
  // hence |delta_phi(r)| should be unimodal
  // use godel section search to find the minimum of |delta_phi(r)| within range [0, max_r]
	double a = 0, b = max_r;
	double k = (sqrt(5.) - 1.) / 2.;
	double xL = b - k * (b - a);
	double xR = a + k * (b - a);
	double EPSILON = 0.001;
	int iter = 0, max_iter = 50;

	while ( fabs(b-a) > EPSILON*(fabs(xR)+fabs(xL)) )
	{
		if ( iter>max_iter ) { cout<<"nothing found within interation!!!"<<endl; r_conv = -9999.; return 0; }

		float dphi_xR, dphi_xL;
		if ( !delta_phi( alpha_e, alpha_p, phi_e, phi_p, xR, z_DC_m_z_ver_e, z_DC_m_z_ver_p, alpha_r_z, dphi_xR ) )
		{ // if delta_phi in xR, [xR, b] will be the new searching range
			a = xR;
	      	xL = b - k * (b - a);
	      	xR = a + k * (b - a);
	      	++iter;
	      	continue;
		}
		if ( !delta_phi( alpha_e, alpha_p, phi_e, phi_p, xL, z_DC_m_z_ver_e, z_DC_m_z_ver_p, alpha_r_z, dphi_xL ) )
		{ // if delta_phi in xL, [xL, b] will be the new searching range, xR will be the new xL
			a = xL;
	      	xL = xR;
	      	xR = a + k * (b - a);
	      	++iter;
	      	continue;
		}
	  	if( TMath::Abs(dphi_xL) < TMath::Abs(dphi_xR) )	
	    { // find minimum of |dphi|
	      	b = xR;
	      	xR = xL;
	      	xL = b - k*(b - a);
	      	++iter;
	    }
	    else
	    {
	      	a = xL;
	      	xL = xR;
	      	xR = a + k * (b - a);
	      	++iter;
	    }
	}

	r_conv = (a + b) / 2.;
	float dphi;
	delta_phi( alpha_e, alpha_p, phi_e, phi_p, r_conv, z_DC_m_z_ver_e, z_DC_m_z_ver_p, alpha_r_z, dphi );
	// cout<<"r_conv "<<r_conv<<" delta_phi "<<dphi<<endl;
	if ( dphi < -9000 ) { r_conv = -9999.; return 0; } // can't do lookup_fit for the whole r range [0, max_r], dphi always -9999., and that's the minimum value

	// if ( r_conv>29 || r_conv<1 )
	// {
	// 	cout<<"******************"<<endl;
	// 	float r_test = 0;
	// 	for (int i = 0; i < 60; ++i)
	// 	{
	// 		delta_phi( alpha_e, alpha_p, phi_e, phi_p, r_test, z_DC_m_z_ver_e, z_DC_m_z_ver_p, alpha_r_z, dphi );
	// 		cout<<"r_test "<<r_test<<" delta_phi "<<dphi<<endl;
	// 		r_test = r_test+0.5;
	// 	}
	// 	cout<<"******************"<<endl;
	// }
	

	return 1;
}


void Reconstruction::findIntersection(MyTrack const* trk1, MyTrack const* trk2, MyPair* pair, float zvertex)
{
	if ( ( trk1->GetAlphaPrime() )*( trk2->GetAlphaPrime() )>0 ) { cout<<"likesign pair!!!"<<endl; return; }

	float r_conv, phi_conv_e, phi_conv_p, theta_conv_e, theta_conv_p;
	if ( trk1->GetAlphaPrime()>0 && trk2->GetAlphaPrime()<0)
	{ 
		find_rconv( trk1->GetAlphaPrime(), trk2->GetAlphaPrime(), trk1->GetPhiDC(), trk2->GetPhiDC(), trk1->GetZDC() - zvertex, trk2->GetZDC() - zvertex, pimpl->alpha_r_z, r_conv );
		project_phi( trk1->GetAlphaPrime(), trk1->GetPhiDC(), r_conv, trk1->GetZDC() - zvertex, pimpl->alpha_r_z, phi_conv_e);
		project_phi( trk2->GetAlphaPrime(), trk2->GetPhiDC(), r_conv, trk2->GetZDC() - zvertex, pimpl->alpha_r_z, phi_conv_p );
		project_theta( trk1->GetAlphaPrime(), r_conv, trk1->GetZDC() - zvertex, pimpl->alpha_r_z, theta_conv_e );
		project_theta( trk2->GetAlphaPrime(), r_conv, trk2->GetZDC() - zvertex, pimpl->alpha_r_z, theta_conv_p );
	}
	else 
	{
		find_rconv( trk2->GetAlphaPrime(), trk1->GetAlphaPrime(), trk2->GetPhiDC(), trk1->GetPhiDC(), trk2->GetZDC() - zvertex, trk1->GetZDC() - zvertex, pimpl->alpha_r_z, r_conv );
		project_phi( trk2->GetAlphaPrime(), trk2->GetPhiDC(), r_conv, trk2->GetZDC() - zvertex, pimpl->alpha_r_z, phi_conv_e);
		project_phi( trk1->GetAlphaPrime(), trk1->GetPhiDC(), r_conv, trk1->GetZDC() - zvertex, pimpl->alpha_r_z, phi_conv_p );
		project_theta( trk2->GetAlphaPrime(), r_conv, trk2->GetZDC() - zvertex, pimpl->alpha_r_z, theta_conv_e );
		project_theta( trk1->GetAlphaPrime(), r_conv, trk1->GetZDC() - zvertex, pimpl->alpha_r_z, theta_conv_p );
	}
	
	pair->SetRPair( r_conv );
	if ( r_conv<-9000. )
	{
		pair->SetPhiElectron( -9999. );
		pair->SetPhiPositron( -9999. );
		pair->SetThetaElectron( -9999. );
		pair->SetThetaPositron( -9999. );
		return;
	}		
	pair->SetPhiElectron( phi_conv_e );
	pair->SetPhiPositron( phi_conv_p );
	pair->SetThetaElectron( theta_conv_e );
	pair->SetThetaPositron( theta_conv_p );
}

TVector3 Reconstruction::findMomentum(MyTrack* trk, float r, float phi_conv, float theta_conv, float zvertex)
{
	TVector3 vec;
	float mom, phi, theta;
	lookup_fit( TMath::Abs(trk->GetAlphaPrime()), r, TMath::Abs(trk->GetZDC() - zvertex), 1, pimpl->alpha_r_z, mom );
	project_phi( trk->GetAlphaPrime(), trk->GetPhiDC(), r, TMath::Abs(trk->GetZDC() - zvertex), pimpl->alpha_r_z, phi );
	project_theta( trk->GetAlphaPrime(), r, trk->GetZDC() - zvertex, pimpl->alpha_r_z, theta );

	float pt = mom*sin(theta);
	vec.SetPtThetaPhi( pt, theta, phi );

	return vec;
}


}
