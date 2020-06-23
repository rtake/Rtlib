# include <vector>
# include <fstream>

# include "Rtlib.h"

using namespace std;

void ConvertXYZtoFRG(int argc, char* argv[]) {
	int natom, i, j, num;
	FILE *frg;
	char buf[256];

	ifstream xyz;
	string infname, outfname, line;
	
	vector< vector<Atom> > mols;
	vector<string> comments;
	vector< vector<int> > fmat;
	vector<int> fvec;

	infname = argv[1];
	outfname = infname + ".frg";

	xyz.open( infname.c_str() );
	while( getline(xyz, line) ) {
		sscanf( line.c_str(), "%d", &natom);
		vector<Atom> m(natom);
		getline( xyz, line ); // cout << line << endl;
		comments.push_back( line );
		// getline(xyz, line);

		for(i = 0;i < natom;i++) {
			Atom a;
			getline(xyz, line);
			a.SetfromString(line);
			m[i] = a;
		}

		fvec = MakeFragment(m, 1.2, 1); // threshold = 1.2
		fmat.push_back(fvec);
		mols.push_back(m);
	}
	xyz.close();

	num = (int)mols.size();
	frg = fopen( outfname.c_str(), "w");
	for(i = 0;i < num;i++) {
		fprintf(frg, "%s\t", comments[i].c_str() ); // cout << comments[i] << endl;
		for(j = 0;j < fmat[i].size();j++) fprintf(frg, "%d\t", fmat[i][j] );
		fprintf(frg, "\n");
	}
	fclose( frg );
}


int main(int argc, char* argv[]) {
	ConvertXYZtoFRG(argc, argv);	
	return 0;
}
