# include "Rtlib.h"

int main( int argc, char* argv[] ) {
	int i, j, num0, num1, natom, **chkmat;
	const int maxsize = 100000;
	FILE *fp;
	char line[4096];
	double threshold = 1.0;
	Atom **mols0, **mols1, a;

	mols0 = new Atom* [maxsize];
	mols1 = new Atom* [maxsize];
	for(i = 0;i < maxsize;i++) { 
		mols0[i] = new Atom [200];
		mols1[i] = new Atom [200];
	}

	fp = fopen( argv[1], "r" );
	num0 = 0;
	while( fgets( line, 4096, fp ) ) {
		sscanf( line, "%d", &natom );
		fgets( line, 4096, fp );

		for(i = 0;i < natom;i++) {
			fgets( line, 4096, fp );
			a.SetfromString( line );
			mols0[num0][i] = a;
		}

		num0++;
	}
	fclose( fp );

	fp = fopen( argv[2], "r" );
	num1 = 0;
	while( fgets( line, 4096, fp ) ) {
		sscanf( line, "%d", &natom );
		fgets( line, 4096, fp );

		for(i = 0;i < natom;i++) {
			fgets( line, 4096, fp );
			a.SetfromString( line );
			mols1[num1][i] = a;
		}

		num1++;
	}
	fclose( fp );

	chkmat = ( int** )malloc( sizeof( int* ) * num0 );
	for(i = 0;i < num0;i++) { chkmat[i] = ( int* )malloc( sizeof( int ) * num1 ); }

	for(i = 0;i < num0;i++) {
		for(j = 0;j < num1;j++) {
			if( DistMat_diff( mols0[i], mols1[j], natom ) < threshold ) { chkmat[i][j] = 1; }
			else { chkmat[i][j] = -1; }
		}
	}

	for(i = 0;i < num0;i++) { for(j = 0;j < num1;j++) { fprintf( stdout, "mol0[%d] - mol1[%d] : %d\n", i, j, chkmat[i][j] ); } }
        for(i = 0;i < num1;i++) { for(j = 0;j < num0;j++) { fprintf( stdout, "mol1[%d] - mol0[%d] : %d\n", i, j, chkmat[j][i] ); } }

	free( chkmat );

	for(i = 0;i < maxsize;i++) {
		delete [] mols0[i];
		delete [] mols1[i];
	}
	delete mols0;
	delete mols1;

	return 0;
}
