# ifndef INCLUDE_GUARD_RTLIB_HPP
# define INCLUDE_GUARD_RTLIB_HPP

# include <iostream>
# include <string>
# include <cstring>
# include <fstream>
# include <stdlib.h>
# include <vector>
# include <cmath>
# include <set>
# include <sstream>
# include <utility>
# include <algorithm>

# include "eigen-3.3.7/Eigen/Dense"
# include "eigen-3.3.7/Eigen/Core"
# include "eigen-3.3.7/Eigen/LU"

using namespace std;
using namespace Eigen;

//////////////////////////////////////////////////////////////////////////////////////////////////

double ang_to_bohr( double d_ang ) { return d_ang / 0.5291772; }
double bohr_to_ang( double d_bohr ) { return d_bohr * 0.5291772; }

class Atom {
	public:
		Atom() : crd(3,0), elm("X") {}; // constructor
		Atom(const Atom &a) : elm(a.elm), crd(a.crd) {} // copy constructor
		vector<double> GetCrd() { return crd; }
		double GetCrd( int i ) { return crd[i]; }
		void SetCrd( int i, double d ) { crd[i] = d; }
		string GetElm() { return elm; }
		void SetElm( string elm ) { this->elm = elm; }
		double radii();
		int mass();
		void angtobohr() { for(int i = 0;i < 3;i++) { crd[i] /= 0.5291772; } }
		void bohrtoang() { for(int i = 0;i < 3;i++) { crd[i] *= 0.5291772; } }

		Atom& operator=( Atom a ) {
			elm = a.elm;
			crd = a.crd;
			return *this;
		}

		int SetfromString(string line) {
			char buf[3]; vector<double> con(3,0);
			int chk = sscanf(line.c_str(),"%s%17lf%17lf%17lf",buf,&crd[0],&crd[1],&crd[2]);

			if(chk == 4) { string name(buf); this->elm = name; return 0; }
			else return -1;
		}

		void Print() { printf("%s%17.12lf%17.12lf%17.12lf\n",elm.c_str(),crd[0],crd[1],crd[2]); }

	private:
		string elm; // element
		vector<double> crd; // coordinate
};


double Atom::radii() {
	if(elm == "H") return 0.32;
	else if(elm == "C") return 0.75;
	else if(elm == "N") return 0.71;
	else if(elm == "O") return 0.63;
	else if(elm == "S") return 1.05;
	else return -99999;
}


int Atom::mass() {
	if(elm == "H") return 1;
	else if(elm == "C") return 12;
	else if(elm == "N") return 14;
	else if(elm == "O") return 16;
	else if(elm == "S") return 32;
	else return -99999;
}


double Dist(Atom a0, Atom a1) {
	double sum = 0;
	vector<double> c0 = a0.GetCrd(), c1 = a1.GetCrd();

	for(int i = 0;i < 3;i++) sum += pow( (c0[i] - c1[i]), 2);
	return sqrt(sum);
}


double Angle( Atom a0, Atom a1, Atom a2 ) {
	int i;
	double innerproduct = 0;
	for(i = 0;i < 3;i++) { innerproduct += ( a0.GetCrd( i ) - a1.GetCrd( i ) ) * ( a2.GetCrd( i ) - a1.GetCrd( i ) ); }
	return acos( innerproduct / ( Dist( a0, a1 ) * Dist( a2, a1 ) ) );
}


double Dihedral( Atom a0, Atom a1, Atom a2, Atom a3 ) {
  Vector3d v0,v1,v2;

  for(int i=0;i<3;i++) { v0(i) = a2.GetCrd(i) - a1.GetCrd(i); }
  for(int i=0;i<3;i++) { v1(i) = a3.GetCrd(i) - a1.GetCrd(i); }
  for(int i=0;i<3;i++) { v2(i) = a0.GetCrd(i) - a1.GetCrd(i); }

  double cos_theta = (v0.cross(v1)).dot(v0.cross(v2))/(v0.cross(v1).norm()*(v0.cross(v2).norm()));
  double sin_theta = v0.dot((v0.cross(v1)).cross(v0.cross(v2)))/(v0.norm()*v0.cross(v1).norm()*(v0.cross(v2).norm()));
  return atan2(sin_theta,cos_theta);
}

//////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////

void ConvertMOLtoDIST( int natom, Atom* m, double* vec ) {
	int index = 0, i;
	for(i = 0;i < natom;i++) { for(int j = i + 1;j < natom;j++, index++) { vec[index] = Dist( m[i], m[j] ); } }
}


int GetPointfromXYZFILE( FILE* fp, Atom** mols_ref, double* enes_ref ) {
	int natom, nref = 0, i;
	char *pt, line[256];

	while( fgets( line, 256, fp ) ) {
		sscanf( line, "%d", &natom );
		fgets( line, 256, fp );
		pt = strstr( line, "," );
		sscanf( pt + 1, "%17lf", &enes_ref[nref] );

		for(i = 0;i < natom && fgets( line, 256, fp );i++) { mols_ref[nref][i].SetfromString( line ); }
		nref++;
	}

	return nref;
}

//////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////

typedef struct InternalCoordinate {
	int natom;
	char **atype; // atom type
	double *r, *a, *d; // radius, angle, dihedral
} InternalCoordinate;


void InternalCoordinate_Malloc( InternalCoordinate *icrd, int num ) {
	icrd = ( InternalCoordinate* )malloc( sizeof( InternalCoordinate ) );
	icrd->atype = ( char** )malloc( sizeof( char* ) * num );
	icrd->r = ( double* )malloc( sizeof( double ) * num );
	icrd->a = ( double* )malloc( sizeof( double ) * num );
	icrd->d = ( double* )malloc( sizeof( double ) * num );
	for(int i = 0;i < num;i++) { icrd->atype[i] = ( char*)malloc( sizeof( char ) * 5 ); }
	icrd->natom = num;
}


void InternalCoordinate_Free( InternalCoordinate *icrd ) {
	if( icrd == NULL ) return;
	for(int i = 0;i < icrd->natom;i++) { free( icrd->atype[i] ); }
	free( icrd->atype  );
	free( icrd->r  );
	free( icrd->a  );
	free( icrd->d  );
	free( icrd );
}

/*
void CartesianToInternal( Atom *mol, InternalCoordinate *icrd ) {
	if( icrd == NULL ) return;
	for(int i = 0;i < icrd->natom;i++) {
		sprintf( icrd->atype[i], "%s", mol[i].GetElm().c_str() );
		icrd->r[i] = icrd->a[i] = icrd->d[i] = -99999;
		if( i > 0 ) icrd->r[i] = Dist( mol[i], mol[0] );
		if( i > 1 ) icrd->a[i] = Angle( mol[i], mol[0], mol[1] );
		if( i > 2 ) icrd->d[i] = Dihedral( mol[i], mol[0], mol[1], mol[2] );
	}
}
*/

void CartesianToInternal( int natom, Atom *mol, double *retv ) {
	const int dfree = 3 * natom - 6;
	if( retv == NULL ) return;
	for(int i = 0;i < dfree;i++) {
		if( i == 0 ) { retv[i] = Dist( mol[1], mol[0] ); }
		else if( i == 1 ) { retv[i] = Dist( mol[2], mol[0] ); }
		else if( i == 2 ) { retv[i] = Angle( mol[2], mol[0], mol[1] ); }
		else if( i %3 == 0 ) { retv[i] = Dist( mol[i/3 + 2], mol[0] ); }
		else if( i %3 == 1 ) { retv[i] = Angle( mol[i/3 + 2], mol[0], mol[1] ); }
		else if( i %3 == 2 ) { retv[i] = Dihedral( mol[i/3 + 2], mol[0], mol[1], mol[2] ); }
	}
}


void InternalToCartesian( int natom, Atom *mol, double *argv ) {
	for(int i = 0;i < natom;i++) {
		for(int j = 0;j < 3;j++) { mol[i].SetCrd( j, 0.0 ); }

		if( i == 1 ) { mol[i].SetCrd( 0, argv[0] ); }
		else if( i == 2 ) {
			mol[i].SetCrd( 0, argv[1] * cos( argv[2] ) );
			mol[i].SetCrd( 1, argv[1] * sin( argv[2] ) );
		} else if( i > 2 ) {
			mol[i].SetCrd( 0, argv[3 * i - 6] * cos( argv[3 * i - 5] ));
			mol[i].SetCrd( 1, argv[3 * i - 6] * sin( argv[3 * i - 5] ) * cos( argv[3 * i - 4]));
			mol[i].SetCrd( 2, argv[3 * i - 6] * sin( argv[3 * i - 5] ) * sin( argv[3 * i - 4]));
		}
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////

void inPESdata(ifstream& ifs, vector< vector<double> >& mat_f, vector<double>& vec_f) {
	// vec_f[i]  \t  mat_f[i][0]  \t  mat_f[i][1]  \t  ...  \n

	string s;
	for(int i = 0;getline(ifs, s);i++) { // cout << s << endl;
		stringstream ssline(s);
		string line;
		vector<double> vals; // vector of distance of each atom pair

		for(int j = 0;getline(ssline,line,'\t');j++) { // cout << "line\t" << line << endl;
			double v;
			int chk = sscanf(line.c_str(),"%17lf",&v);
			if(chk <= 0) printf("sscanf in inPESdata() failed, i : %d, j : %d\n",i,j);

			if(j == 0) vec_f[i] = v; // vec_f.push_back(v);			
			else if(j > 0) vals.push_back(v);
		}
	
		mat_f[i] = vals; // mat_f.push_back(vals);
	}
	// printf("inPESdata() ok\n");
}


void outPESdata(ofstream& ofs, vector< vector<double> > mat_f, vector<double> vec_f) {
	// vec_f[i]  \t  mat_f[i][0]  \t  mat_f[i][1]  \t  ...  \n

	int n = (int)mat_f.size(); // vector of (vector of distance)
	
	if( n != (int)vec_f.size() ) { cout << "number of rows not matched !\n"; return; }
	// else printf("natom chk --> ok (n : %d)\n",n);

	for(int i = 0;i < n;i++) {
		ofs << vec_f[i]; cout << vec_f[i];
		for(int j = 0;j < (int)mat_f[i].size();j++) { ofs << "\t" << mat_f[i][j]; cout << "\t" << mat_f[i][j]; }
		ofs << endl; cout << endl;
	}

	printf("outPESdata() ok\n");
}

//////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////

int BondJudge(Atom a0, Atom a1, double threshold = 1.5) {
	if( Dist(a0, a1) < ( a0.radii() + a1.radii() ) * threshold ) return 1;
	else return -1;
}

/*
vector<int> MakeFragment(vector<Atom> mol, double threshold = 1.5) {
  int i ,j, k;
  const int natom = (int)mol.size();
  vector<int> frg, fvec(natom, 0), f(natom), root(natom); // root
  vector< vector<int> > mat(natom, vector<int>(natom, -1) );

  // 1. Make distance matrix

  for(i = 0;i < natom;i++) {
    for(j = i + 1;j < natom;j++) {
      if( BondJudge(mol[i], mol[j], threshold) == 1 ) {
        mat[i][j] = 1;
        mat[j][i] = 1;
      }
    }
  }
	
  for(i = 0;i < natom;i++) { f[i] = i; root[i] = i; }

  for(i = 0;i < natom;i++) {
    for(j = i + 1;j < natom;j++) {
      if( mat[i][j] == 1 ) {
        for(k = 0;k < natom;k++) {
          if( k == j ) continue;
          if( root[k] == root[j] ) { root[k] = root[i]; }
        }
        root[ root[j] ] = root[i];
        root[j] = root[i];
      }
    }
  }

  // 2. Calc. Mass of each fragments

  for(i = 0;i < natom;i++) { fvec[ root[i] ] += mol[i].mass(); }
  for(i = 0;i < natom;i++) { if( fvec[i] > 0 ) frg.push_back( fvec[i] ); }
  sort( frg.begin(), frg.end(), greater<int>() );

  for(i = 0;i < natom;i++) { fprintf( stdout, "root[%d] : %d\t", i, root[i] ); }
  fprintf( stdout, "\n" );

  return frg;
}
*/

vector<int> MakeFragment(vector<Atom> mol, double threshold = 1.5, int type=0) {
  int natom = (int)mol.size();
  vector<int> root(natom), frg, fvec(natom, 0), f(natom);
  vector< vector<int> > dmat(natom,vector<int>(natom,-1));

  for(int i=0;i<natom;i++) {
    for(int j=i+1;j<natom;j++) {
      if(BondJudge(mol[i], mol[j], threshold) == 1) {
        dmat[i][j] = 1;
        dmat[j][i] = 1;
      }
    }
  }

  for(int i=0;i<natom;i++) { root[i] = i; }

  for(int i=0;i<natom;i++) {
    for(int j=i+1;j<natom;j++) {
      if(dmat[i][j] != 1) { continue; } // skip no bond

      for(int k=0;k<natom;k++) {
        if(k==j) { continue; } // same atom
        if(root[k] == root[j]) { root[k] = root[i]; }
      }
      root[ root[j] ] = root[i];
      root[j] = root[i];
    }
  }

  if(type==0) { return root; }


  for(int i=0;i<natom;i++) { fvec[ root[i] ] += mol[i].mass(); }
  for(int i=0;i<natom;i++) { if( fvec[i] > 0 ) frg.push_back( fvec[i] ); }
  sort( frg.begin(), frg.end(), greater<int>() );

  for(int i=0;i<natom;i++) { fprintf( stdout, "root[%d] : %d\t", i, root[i] ); }
  fprintf( stdout, "\n" );

  return frg;
}

//////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////

int ConvertLISTLOGtoXYZ(int argc, char* argv[]) {
	int i, j, num = 0, nmode;
	FILE *fp;
	double ene, spin, zpve, evalue;
	
	string comment, line, infname, outfname;
	ifstream list;
	vector<double> evvec;
	vector<string> svec;
	pair<int,int> connection;

	list.open( argv[1] );
	infname = argv[1];
	outfname = infname + ".xyz";
	fp = fopen( outfname.c_str(), "w");

	while( getline( list, line ) ) { if( strstr( line.c_str(), "# Geometry") ) break; }

	do {
		// Load start

		svec.clear();
		evvec.clear();
		ene = spin = zpve = evalue = -99999;
		connection.first = -1; connection.second = -1;

		comment = line;

		while( getline( list, line ) ) {
			if( strstr(line.c_str(), "Energy") ) { break; }
			svec.push_back( line );
		}
		
		if( strstr( line.c_str(), "Energy") ) { sscanf( line.c_str() + 11, "%17lf", &ene); getline( list, line ); }
		if( strstr( line.c_str(), "Spin") ) { sscanf( line.c_str() + 11, "%17lf", &spin); getline( list, line ); }
		if( strstr( line.c_str(), "ZPVE") ) { sscanf( line.c_str() + 11, "%17lf", &zpve); getline( list, line ); }
		if( strstr( line.c_str(), "eigenvalues") ) {
			sscanf( line.c_str() + 34, "%d", &nmode); getline( list, line );
			evvec.resize( nmode );
			for(i = 0;i < nmode; getline( list, line ) ) {
				for(j = 0;j < 5 && i < nmode;j++, i++) sscanf( line.c_str() + 1 + 14*j, "%11lf", &evvec[i] );
			}
		}
		if( strstr( line.c_str(), "CONNECTION") ) { sscanf( line.c_str() + 13, "%d - %d", &connection.first, &connection.second ); }

		// Load end

		// Print start

		fprintf(fp, "%d\n", (int)svec.size() );
		fprintf(fp, "%s,Energy=%17.12lf,Spin(**2)=%17.12lf,ZPVE=%17.12lf,", comment.c_str(), ene, spin, zpve);
//		fprintf(fp, "nmode=%d/", nmode);
//		for(i = 0;i < nmode;i++) { fprintf(fp, "%10.9lf/", evvec[i] ); }
//		fprintf(fp, "CONNECTION:%d-%d/\n", connection.first, connection.second );
		fprintf(fp,"\n");
		for(i = 0;i < (int)svec.size();i++) { fprintf(fp, "%s\n", svec[i].c_str() ); }

		// Print end

		num++;
	} while( getline( list, line ) );

	list.close();

	return num;
}


int ConvertIRCLOGtoXYZ(int argc, char* argv[]) {
	char type[256], line[256], xyz[256];
	double ene, spn; // energy, spin
	int i, j, natom = 0, chk = 0; // number of step
	FILE *fp;

	vector< vector<Atom> > all_mols, irc_mols;
	vector< double > v_ene, v_spn;
	vector< vector<Atom> >::iterator itr_m;
	vector< double >::iterator itr_e, itr_s;
	vector< Atom > fwd, bck;

	sprintf( xyz, "%s.xyz", argv[1]);

	fp = fopen( argv[1], "r" );
	while( fgets( line, 256, fp ) ) {
		sprintf( type, "OFF" );
		if( strstr( line, "INITIAL STRUCTURE") ) { sprintf( type, "TS" ); }
		else if( strstr( line, "# STEP") ) { sprintf( type, "IRC" ); }
		else if( strstr( line, "IRC FOLLOWING (FORWARD)") ) { chk = 1; } // chk == 1 --> forward
		else if( strstr( line, "IRC FOLLOWING (BACKWARD)") ) { chk = -1; } // chk == -1 --> backward
		else if( strstr( line, "Optimized") ) { sprintf( type, "OPT" ); }

		if( strstr( type, "OFF" ) ) continue;

		vector< Atom > m;
		Atom a;

		while( fgets( line, 256, fp ) ) {
			if( a.SetfromString( line ) != 0 ) break;
			m.push_back(a);
		}
		sscanf( line + 11, "%17lf", &ene );
		fgets( line, 256, fp );
		sscanf( line + 11, "%17lf", &spn );

		if( chk <= 0 ) { // TS, BACKWARD
			all_mols.push_back( m );
			v_ene.push_back( ene );
			v_spn.push_back( spn );
			itr_m = all_mols.begin();
			itr_e = v_ene.begin();
			itr_s = v_spn.begin();
		} else { // FORWARD
			itr_m = all_mols.insert( itr_m, m );
			itr_e = v_ene.insert( itr_e, ene );
			itr_s = v_spn.insert( itr_s, spn );
		}

		if( strstr( type, "OPT") ) {
			if( chk > 0 ) { fwd = m; }
			else if( chk < 0 ) { bck = m; }
		}


	}
	fclose( fp );

	if( fwd.size() == all_mols[0].size() ) { all_mols[0] = fwd; }
	else if( bck.size() == all_mols[0].size() ) { all_mols[ (int)all_mols.size() - 1 ] = bck; }

	fp = fopen( xyz, "w" );
	natom = (int)all_mols[0].size();
	for(i = 0;i < (int)all_mols.size();i++) {
		fprintf( fp, "%d\nStr. %d/%17.12lf/%17.12lf\n", natom, i, v_ene[i], v_spn[i] );
		for(j = 0;j < natom;j++) {
			fprintf( fp, "%s\t%17.12lf\t%17.12lf\t%17.12lf\n",
			all_mols[i][j].GetElm().c_str(), all_mols[i][j].GetCrd()[0], all_mols[i][j].GetCrd()[1], all_mols[i][j].GetCrd()[2] );
		}
	}
	fclose( fp );

	return 0;
}


int ConvertXYZtoDIST(int argc, char* argv[]) {
	FILE *fp;
	char *pt, line[256], file[256];
	int natom, i, j, k;
	double e;

	vector< vector< Atom > > mols;
	vector< double > vec_e;
	Atom a;

	sprintf( file, "%s.dist", argv[1] );
	fp = fopen( argv[1], "r");
	while( fgets( line, 256, fp ) ) {
		sscanf( line, "%d", &natom);
		fgets( line, 256, fp );
		pt = strstr( line, "/" ); 
		sscanf( pt + 1, "%lf", &e );

		vector< Atom > m( natom );
		for(i = 0;i < natom;i++) {
			fgets( line, 256, fp );
			a.SetfromString( line );
			m[i] = a;
		}
		mols.push_back( m );
		vec_e.push_back( e );
	}
	fclose( fp );

	fp = fopen( file, "w" );
	for(i = 0;i < (int)mols.size();i++) {
		fprintf( fp, "%17.12lf\t", vec_e[i] );
		for(j = 0;j < natom;j++) {		
			for(k = j + 1;k < natom;k++) {
				fprintf( fp, "%17.12lf\t", Dist( mols[i][j], mols[i][k] ) );
			}
		}
		fprintf( fp, "\n" );
	}
	fclose( fp );

	return 0;
}


int ConvertLUPOUTttoXYZ(int argc, char* argv[]) {
	FILE *fpxyz;
	double ene, spn;
	int i;

	ifstream log;
	string infname, outfname, line, comment;
	vector< string > svec;
	
	infname = argv[1];
	outfname = infname + ".xyz";
	fpxyz = fopen( outfname.c_str(), "w" );
	log.open( argv[1] );
	while( getline( log, line ) ) {
		if( !strstr( line.c_str(), "# NODE") ) { continue; }

		svec.clear();
		comment = line;
		while( getline( log, line ) ) {
			if( strstr( line.c_str(), "Threshold" ) ) { break; }
			svec.push_back( line );
		}

		while( getline( log, line ) ) {
			if( strstr( line.c_str(), "ENERGY" ) ) {
				sscanf( line.c_str() + 22, "%17lf", &ene);
			} else if( strstr( line.c_str(), "Spin(**2)" ) ) {
				sscanf( line.c_str() + 22, "%17lf", &spn);
				break;
			}
		}

		fprintf( fpxyz, "%d\n", (int)svec.size() );
		fprintf( fpxyz, "%s/%17.12lf/%17.12lf\n", comment.c_str(), ene, spn );
		for(i = 0;i < (int)svec.size();i++) { fprintf(fpxyz, "%s\n", svec[i].c_str() ); }

	}
	log.close();
	fclose( fpxyz );

	return 0;
}


int ConvertMINtoXYZ(int argc, char* argv[]) {
	FILE *fpxyz;
	double ene;
	int i;
	const char *pt;

	ifstream log;
	string infname, outfname, line, comment;
	vector< string > svec;

	infname = argv[1];
	outfname = infname + ".xyz";
	fpxyz = fopen( outfname.c_str(), "w" );
	log.open( argv[1] );
	while( getline( log, line ) ) {
		if( !strstr( line.c_str(), "# ITR") ) { continue; }

		svec.clear();
		comment = line;
		while( getline( log, line ) ) {
			if( strstr( line.c_str(), "Item" ) ) { break; }
			svec.push_back( line );
		}

		while( getline( log, line ) ) {
			if( strstr( line.c_str(), "ENERGY" ) ) {
				pt = strstr( line.c_str(), "ENERGY" );
				sscanf( pt + 10, "%17lf", &ene);
				break;
			}
		}

		fprintf( fpxyz, "%d\n", (int)svec.size() );
		fprintf( fpxyz, "%s/%17.12lf\n", comment.c_str(), ene );
		for(i = 0;i < (int)svec.size();i++) { fprintf(fpxyz, "%s\n", svec[i].c_str() ); }

	}
	log.close();
	fclose( fpxyz );

	return 0;
}

//////////////////////////////////////////////////////////////////////////////////////////////////

int Combination(int n, int r) {
	if(n == r) return 1;
	else if(r == 0) return 1;
	else if(r == 1) return n;
	else return Combination(n - 1, r - 1) + Combination(n - 1, r);

}


void MakeCombination( int size, int order, int **mat ) {
	int i, j, k, index, **mat0, size0;

	index = 0;
	if( size <= 1 ) {
		for(i = 0;i < order + 1;i++) { mat[i][index] = i; }
	} else {
		for(i = 0;i < order + 1;i++) {
			size0 = Combination( size - 1 + i, i ); // nrow
			mat0 = ( int** )malloc( sizeof( int* ) * ( size0 ) ); // row
			for(j = 0;j < size0;j++) { mat0[j] = ( int* )malloc( sizeof( int* ) * ( size - 1 ) ); } // line
			MakeCombination( size - 1, i, mat0 );

			for(j = 0;j < size0;j++, index++) {
				mat[index][0] = order - i;
				for(k = 0;k < size - 1;k++) { mat[index][k + 1] = mat0[j][k]; }
			}
			
			for(j = 0;j < size0;j++) { free( mat0[j] ); }
			free( mat0 );
		}
	}

}

//////////////////////////////////////////////////////////////////////////////////////////////////

double DistMat_diff( Atom *mol0, Atom *mol1, int natom ) {
	int i, j;
	double diff = 0;
	for(i = 0;i < natom;i++) { for(j = i + 1;j < natom;j++) { diff += abs( Dist( mol0[i],mol0[j] ) - Dist( mol1[i],mol1[j] ) ); } }
	return diff;
}

//////////////////////////////////////////////////////////////////////////////////////////////////

# endif
