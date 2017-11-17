
#include "DnbsClassGamma.h"
#include "LyonLib.h"

DnbsClassGamma::DnbsClassGamma(void)
{
}

DnbsClassGamma::~DnbsClassGamma(void)
{
}

bool DnbsClassGamma::genHaarBases_inProd(double minInnProd) {
	// Divide whole dictionary to clusters;
	// Each cluster has one central base vector;
	// The inner product between any of the features 
	// and the central vector is larger than 'minInnProd'.

	// Initialization

	inprod_limit = minInnProd;

	// Clear previous data
	indexOfCentBases.clear();
	for (int i = 0; i < (int)memberOfGroup.size(); ++ i)
		memberOfGroup[i].clear();
	memberOfGroup.clear();

	// Load if saved
	char filename[256];
	sprintf(filename, "w%dh%dprod%.6lf.dnbsdic", imagW, imagH, minInnProd);
	FILE * fin = fopen(filename, "r");
	int numGroups;
	if (fin != NULL) {
		fscanf(fin, "%d", &numGroups);
		printf("# Clusters = %d\n", numGroups);
		for (int i = 0; i < numGroups; ++ i) {
			int cent, numMembers;
			fscanf(fin, "%d%d", &cent, &numMembers);
			indexOfCentBases.push_back(cent);
			memberOfGroup.push_back(vector<int>());
			for (int j = 0; j < numMembers; ++ j) {
				int curMember;
				fscanf(fin, "%d", &curMember);
				memberOfGroup[i].push_back(curMember);
			}
		}
		fclose(fin);
		return true;
	}

	// Compute if not saved

	int N = listOfBases.size();
	int remainN = N;
	int * locIndex = new int[N];
	int * remainIndex = new int[N];
	for (int i = 0; i < N; ++ i)
		locIndex[i] = remainIndex[i] = i;

	while (remainN > 0) {
		int ind = (int)(rand() / (double)RAND_MAX * remainN);
#ifdef __DEBUG
		if (indexOfCentBases.size() < 10) {
			printf("remainN = %d\n", remainN);
			printf("ind = %d\n", ind);
		}
#endif
		indexOfCentBases.push_back(remainIndex[ind]);
		memberOfGroup.push_back(vector<int>());
		vector<int> & locMemOfGroup = memberOfGroup[memberOfGroup.size() - 1];

		HaarBase & locHaarBase = listOfBases[remainIndex[ind]];
		int oW = locHaarBase.w;
		int oH = locHaarBase.h;
		// clear up nearby features
		swap(remainIndex[ind], remainIndex[remainN - 1]);
		swap(locIndex[remainIndex[ind]], locIndex[remainIndex[remainN - 1]]);
		-- remainN;

		for (int i = 0; i < remainN; ) {
			int ii = remainIndex[i];
			if (locHaarBase.innerProd(listOfBases[ii]) >= minInnProd) {
				locMemOfGroup.push_back(ii);
				swap(remainIndex[i], remainIndex[remainN - 1]);
				swap(locIndex[remainIndex[i]], locIndex[remainIndex[remainN - 1]]);
				-- remainN;
				continue;
			}
			++ i;
		}
	}

	printf("Hierarchy generated: %d clusters.\n", indexOfCentBases.size());

	delete [] locIndex;
	delete [] remainIndex;

	// Save results
	FILE * fout = fopen(filename, "w+");
	fprintf(fout, "%d\n", indexOfCentBases.size());
	for (int i = 0; i < (int)indexOfCentBases.size(); ++ i) {
		fprintf(fout, "%d %d", indexOfCentBases[i], memberOfGroup[i].size());
		for (int j = 0; j < (int)memberOfGroup[i].size(); ++ j) {
			fprintf(fout, " %d", memberOfGroup[i][j]);
		}
		fprintf(fout, "\n");
	}
	fclose(fout);

	return true;
}

void DnbsClassGamma::computeDNBS_v1_branchbound(
	int MaxNumBases,
	double lambda,
	vector<CvMat*> foreground,
	vector<CvMat*> background,
	double inprod_limit) {

		// Supposing frame is gray-level image with 3 channels

		// Supposing 'GenHaarBases_innerprod' completed.

		// Record the width and height
		int width = foreground[0]->width;
		int height = foreground[0]->height;

		int num_foregrounds = foreground.size();
		int num_backgrounds = background.size();

		// Definitions
		// K is the number of selected bases
		int K;
		// N is the total number of candidate bases
		int N = listOfBases.size();

		// B and D are two arrays of parameters in deciding the
		// optimized index of the next to-be-selected binary base.
		// i.e. selecting the one that maximising E_n(=|BF_n||BG_n|/D_n)
		vector<NBS_MAT_TYPE*> BF;
		for (int i = 0; i < num_foregrounds; ++ i) {
			NBS_MAT_TYPE * tempBF
				= new NBS_MAT_TYPE[N];
			BF.push_back(tempBF);
		}

		vector<NBS_MAT_TYPE*> BG;
		for (int i = 0; i < num_backgrounds; ++ i) {
			NBS_MAT_TYPE * tempBG
				= new NBS_MAT_TYPE[N];
			BG.push_back(tempBG);
		}
		NBS_MAT_TYPE * D = new NBS_MAT_TYPE[N];

		// Indicating whether the base is selected to form the subspace
		bool * selected = new bool[N];

		// Index of each selected binary bases
		vector<int> L;
		// Projection Coefficients
		vector<NBS_MAT_TYPE> Coeff;
		vector< vector<double> > Coeff_fg, Coeff_bg;

		// Phi's and their integral image
		vector<CvMat*> Phi;
		vector<CvMat*> intPhi;

		// Inner product of each Phi[K] and image F(i.e. matPatch)
		vector< vector<NBS_MAT_TYPE> > innerPhiF;
		for (int i = 0; i < num_foregrounds; ++ i) 
			innerPhiF.push_back(vector<NBS_MAT_TYPE>());
		// Inner product of each Phi[K] and image G
		vector< vector<NBS_MAT_TYPE> > innerPhiG;
		for (int i = 0; i < num_backgrounds; ++ i)
			innerPhiG.push_back(vector<NBS_MAT_TYPE>());

		// Beta's and their integral image
		vector<CvMat*> Beta;
		vector<CvMat*> intBeta;

		vector<double> fres, bres;

		// Initialization
		// Reset memories
		memset(selected, 0, sizeof(bool) * N);
		// Computing the integral image of the original image
		vector<CvMat*> integralF;
		for (int i = 0; i < num_foregrounds; ++ i)
			integralF.push_back(calcIntegral(foreground[i]));
		vector<CvMat*> integralG;
		for (int i = 0; i < num_backgrounds; ++ i)
			integralG.push_back(calcIntegral(background[i]));

		/// Representative features from the dictionary
		int sampN = indexOfCentBases.size();

		int * progress = new int[N];

		// Process of initializing parameters and get the first choice
		int tempL = 0;
		double maxFunc = -inf;
		double curFunc;
		for(int i = 0; i < N; ++ i) {
			progress[i] = 0;
			for (int j = 0; j < num_foregrounds; ++ j)
				BF[j][i] = listOfBases[i].innerProd(integralF[j]);
			for (int j = 0; j < num_backgrounds; ++ j)
				BG[j][i] = listOfBases[i].innerProd(integralG[j]);
			D[i] = 1;

			curFunc = 0;
			for (int j = 0; j < num_foregrounds; ++ j)
				curFunc += BF[j][i] * BF[j][i] / num_foregrounds;
			for (int j = 0; j < num_backgrounds; ++ j)
				curFunc -= BG[j][i] * BG[j][i] / num_backgrounds * lambda;

			if (curFunc > maxFunc) {
				tempL = i;
				maxFunc = curFunc;
			}
		}
		L.push_back(tempL);
		selected[tempL] = true;

		// Computing parameters
		Phi.push_back(listOfBases[L[0]].toMatrix(height, width));
		intPhi.push_back(calcIntegral(Phi[0]));
		Beta.push_back(cvCloneMat(Phi[0]));
		intBeta.push_back(calcIntegral(Beta[0]));
		for (int i = 0; i < num_foregrounds; ++ i) {
			Coeff_fg.push_back(vector<double>());
			Coeff_fg[i].push_back(BF[i][L[0]]);
		}
		for (int i = 0; i < num_backgrounds; ++ i) {
			Coeff_bg.push_back(vector<double>());
			Coeff_bg[i].push_back(BG[i][L[0]]);
		}

		//Coeff.push_back(BF[L[0]]);
		for (int i = 0; i < num_foregrounds; ++ i)
			innerPhiF[i].push_back(innerProd(Phi[0], foreground[i]));
		for (int i = 0; i < num_backgrounds; ++ i)
			innerPhiG[i].push_back(innerProd(Phi[0], background[i]));

		for (int i = 0; i < num_foregrounds; ++ i) {
			fres.push_back(0);
			for (int r = 0; r < foreground[i]->rows; ++ r)
				for (int c = 0; c < foreground[i]->cols; ++ c) {
					double curr = CV_MAT_ELEM(*foreground[i], double, r, c)
						- Coeff_fg[i][0] * listOfBases[L[0]].element(r, c, imagW);
					fres[i] = fres[i] + curr * curr;
				}
				fres[i] = sqrt(fres[i]);
		}

		for (int i = 0; i < num_backgrounds; ++ i) {
			bres.push_back(0);
			for (int r = 0; r < background[i]->rows; ++ r)
				for (int c = 0; c < background[i]->cols; ++ c) {
					double curr = CV_MAT_ELEM(*background[i], double, r, c)
						- Coeff_bg[i][0] * listOfBases[L[0]].element(r, c, imagW);
					bres[i] = bres[i] + curr * curr;
				}
				bres[i] = sqrt(bres[i]);
		}

		int curtime = clock();

		double * sampScore = new double[sampN];
		double * bounds = new double[sampN];

		double vardelta = sqrt(2. - 2. * inprod_limit);

		K = 1;
		// Total Time Complexity is O(KN + K^2WH) = O(KWH(WH+K))
		while (K < MaxNumBases) {
			/**
			xmlDocPtr doc = xmlNewDoc(BAD_CAST"1.0");
			xmlNodePtr root_node = xmlNewNode(NULL, BAD_CAST"staterror");
			xmlDocSetRootElement(doc, root_node);
			static char buf[255];*/

			// Selecting the best fit Base
			maxFunc = -inf;
			tempL = -1;
			double tempPhiBase;
			int optIndex;

			for(int i = 0; i < sampN; ++ i) {
				//if (!selected[i]) {
				int ii = indexOfCentBases[i];
				progress[ii] = K;
				tempPhiBase = listOfBases[ii].innerProd(intPhi[K - 1]);
				for (int j = 0; j < num_foregrounds; ++ j)
					BF[j][ii] -= innerPhiF[j][K - 1] * tempPhiBase;
				for (int j = 0; j < num_backgrounds; ++ j)
					BG[j][ii] -= innerPhiG[j][K - 1] * tempPhiBase;
				D[ii] -= tempPhiBase * tempPhiBase;

				curFunc = 0;
				for (int j = 0; j < num_foregrounds; ++ j)
					curFunc += BF[j][ii] * BF[j][ii] / num_foregrounds;
				for (int j = 0; j < num_backgrounds; ++ j)
					curFunc -= BG[j][ii] * BG[j][ii] / num_backgrounds * lambda;

				curFunc = fabs(curFunc) < sqeps ? 0 : (curFunc / D[ii]);

				sampScore[i] = curFunc;

				/**
				xmlNodePtr node = xmlNewNode(NULL, BAD_CAST "base");
				xmlAddChild(root_node, node);
				xmlNewProp(node, BAD_CAST "index", BAD_CAST itoa(i, buf, 255));
				xmlNewProp(node, BAD_CAST "func", BAD_CAST double2string(curFunc, buf, 255)); */

				if (curFunc > maxFunc) {
					maxFunc = curFunc;
					tempL = ii;
					optIndex = i;
				}

				// Calculate Bounds
				double denom = D[ii] - 2 * vardelta * sqrt(D[ii]);
				if (denom < eps) {
					bounds[i] = inf;
				} else {
					double left_term = 0;
					for (int j = 0; j < num_foregrounds; ++ j)
						left_term += lyon::sqr(BF[j][ii] + vardelta * fres[j]);
					double right_term = 0;
					for (int j = 0; j < num_backgrounds; ++ j)
						right_term += lyon::sqr(MAX(BG[j][ii] - vardelta * bres[j], 0));
					bounds[i] = left_term / num_foregrounds / denom
						- lambda * right_term / num_backgrounds / denom;
					printf("curFunc = %.6lf, bound = %.6lf\n", curFunc, bounds[i]);
				}
			}

			/**
			/// Save XML documents
			//Dumping document to stdio or file
			xmlSaveFormatFileEnc((string("exp/") + itoa(K, buf, 255) + string(".xml")).c_str(), doc, "UTF-8", 1);
			xmlFreeDoc(doc);
			xmlCleanupParser();
			xmlMemoryDump();//debug memory for regression tests */

			if (tempL == -1) break;

			int num_of_hits = 0;

			for (int ind = 0; ind < sampN; ++ ind)
				if (bounds[ind] > maxFunc) {

					++ num_of_hits;

					for (int i = 0; i < (int)memberOfGroup[ind].size(); ++ i) {
						int ii = memberOfGroup[ind][i];

						for (int k = progress[ii]; k < K; ++ k) {
							tempPhiBase = listOfBases[ii].innerProd(intPhi[k]);
							for (int j = 0; j < num_foregrounds; ++ j)
								BF[j][ii] -= innerPhiF[j][k] * tempPhiBase;
							for (int j = 0; j < num_backgrounds; ++ j)
								BG[j][ii] -= innerPhiG[j][k] * tempPhiBase;
							D[ii] -= tempPhiBase * tempPhiBase;
						}

						progress[ii] = K;

						curFunc = 0;
						for (int j = 0; j < num_foregrounds; ++ j)
							curFunc += BF[j][ii] * BF[j][ii] / num_foregrounds;
						for (int j = 0; j < num_backgrounds; ++ j)
							curFunc -= BG[j][ii] * BG[j][ii] / num_backgrounds * lambda;

						curFunc = fabs(curFunc) < sqeps ? 0 : (curFunc / D[ii]);

						if (curFunc > maxFunc) {
							maxFunc = curFunc;
							tempL = ii;
						}

					}
				}

				printf("number of hits = %d\n", num_of_hits);

				L.push_back(tempL);
				selected[tempL] = true;

				// Computing parameters
				CvMat * tempMat = calcPhi(L[K], Phi, intPhi);       // O(KWH)
				Phi.push_back(calcMatDiv(tempMat, sqrt(D[L[K]])));  // O(WH)
				Beta.push_back(calcMatDiv(tempMat, D[L[K]]));       // O(WH)
				intBeta.push_back(calcIntegral(Beta[K]));           // O(WH)
				cvReleaseMat(&tempMat);
				intPhi.push_back(calcIntegral(Phi[K]));             // O(WH)
				for (int i = 0; i < num_foregrounds; ++ i)
					innerPhiF[i].push_back(innerProd(Phi[K], foreground[i]));
				for (int i = 0; i < num_backgrounds; ++ i)
					innerPhiG[i].push_back(innerProd(Phi[K], background[i]));

				// Post-processing at the end of each selection
				for (int i = 0; i < num_foregrounds; ++ i)
					Coeff_fg[i].push_back(BF[i][L[K]] / D[L[K]]);
				for (int i = 0; i < num_backgrounds; ++ i)
					Coeff_bg[i].push_back(BG[i][L[K]] / D[L[K]]);
				NBS_MAT_TYPE temp;
				for(int i = 0; i < K; ++ i) {   // O(KWH)
					temp = listOfBases[L[K]].innerProd(intBeta[i]);
					for(int r = 0; r < Beta[i]->rows; ++ r)
						for(int c = 0; c < Beta[i]->cols; ++ c)
							CV_MAT_ELEM(*Beta[i], NBS_MAT_TYPE, r, c)
							-= CV_MAT_ELEM(*Beta[K], NBS_MAT_TYPE, r, c) * temp;
					// 'temp' is real, so its conjugate is itself;
					for (int j = 0; j < num_foregrounds; ++ j)
						Coeff_fg[j][i] -= temp * Coeff_fg[j][K];
					for (int j = 0; j < num_backgrounds; ++ j)
						Coeff_bg[j][i] -= temp * Coeff_bg[j][K];
					cvReleaseMat(&intBeta[i]);
					intBeta[i] = calcIntegral(Beta[i]);
				}

				for (int i = 0; i < num_foregrounds; ++ i) {
					fres[i] = 0;
					for (int r = 0; r < foreground[i]->rows; ++ r)
						for (int c = 0; c < foreground[i]->cols; ++ c) {
							double curr = CV_MAT_ELEM(*foreground[i], double, r, c);
							for (int j = 0; j < K; ++ j) {
								curr -= Coeff_fg[i][j] * listOfBases[L[j]].element(r, c, imagW);
							}
							fres[i] = fres[i] + curr * curr;
						}
						fres[i] = sqrt(fres[i]);
				}

				for (int i = 0; i < num_backgrounds; ++ i) {
					bres[i] = 0;
					for (int r = 0; r < background[i]->rows; ++ r)
						for (int c = 0; c < background[i]->cols; ++ c) {
							double curr = CV_MAT_ELEM(*background[i], double, r, c);
							for (int j = 0; j < K; ++ j) {
								curr -= Coeff_bg[i][j] * listOfBases[L[j]].element(r, c, imagW);
							}
							bres[i] = bres[i] + curr * curr;
						}
						bres[i] = sqrt(bres[i]);
				}

				// Next iteration
				++ K;
		}

		lyon::settime((clock() - curtime) / (double)CLOCKS_PER_SEC);

		// Copy Answers
		this->indexOfBases.clear();
		for(int i = 0; i < (int)L.size(); ++ i)
			this->indexOfBases.push_back(listOfBases[L[i]]);
		CvMat * intTarget = calcIntegral(target);
		this->coefMatrix();
		this->coeffOfBases = compCoeff(intTarget);

		//demo();

		// Release Memories
		delete [] progress;
		delete [] D;
		delete [] selected;
		for (int i = 0; i < num_foregrounds; ++ i) {
			delete [] BF[i];
		}
		for (int i = 0; i < num_backgrounds; ++ i) {
			delete [] BG[i];
		}
		for(int i = 0; i < (int)Phi.size(); ++ i) {
			cvReleaseMat(&Phi[i]);
			cvReleaseMat(&intPhi[i]);
			cvReleaseMat(&Beta[i]);
			cvReleaseMat(&intBeta[i]);
		}
}

void DnbsClassGamma::computeDNBS_v1_cluster(
	int MaxNumBases,
	double lambda,
	vector<CvMat*> foreground,
	vector<CvMat*> background,
	double ratio) {

		// Supposing frame is gray-level image with 3 channels

		// Supposing 'GenHaarBases_innerprod' completed.

		// Record the width and height
		int width = foreground[0]->width;
		int height = foreground[0]->height;

		int num_foregrounds = foreground.size();
		int num_backgrounds = background.size();

		// Definitions
		// K is the number of selected bases
		int K;
		// N is the total number of candidate bases
		int N = listOfBases.size();

		// B and D are two arrays of parameters in deciding the
		// optimized index of the next to-be-selected binary base.
		// i.e. selecting the one that maximising E_n(=|BF_n||BG_n|/D_n)
		vector<NBS_MAT_TYPE*> BF;
		for (int i = 0; i < num_foregrounds; ++ i) {
			NBS_MAT_TYPE * tempBF
				= new NBS_MAT_TYPE[N];
			BF.push_back(tempBF);
		}

		vector<NBS_MAT_TYPE*> BG;
		for (int i = 0; i < num_backgrounds; ++ i) {
			NBS_MAT_TYPE * tempBG
				= new NBS_MAT_TYPE[N];
			BG.push_back(tempBG);
		}
		NBS_MAT_TYPE * D = new NBS_MAT_TYPE[N];

		// Indicating whether the base is selected to form the subspace
		bool * selected = new bool[N];

		// Index of each selected binary bases
		vector<int> L;
		// Projection Coefficients
		vector<NBS_MAT_TYPE> Coeff;
		vector< vector<double> > Coeff_fg, Coeff_bg;

		// Phi's and their integral image
		vector<CvMat*> Phi;
		vector<CvMat*> intPhi;

		// Inner product of each Phi[K] and image F(i.e. matPatch)
		vector< vector<NBS_MAT_TYPE> > innerPhiF;
		for (int i = 0; i < num_foregrounds; ++ i) 
			innerPhiF.push_back(vector<NBS_MAT_TYPE>());
		// Inner product of each Phi[K] and image G
		vector< vector<NBS_MAT_TYPE> > innerPhiG;
		for (int i = 0; i < num_backgrounds; ++ i)
			innerPhiG.push_back(vector<NBS_MAT_TYPE>());

		// Beta's and their integral image
		vector<CvMat*> Beta;
		vector<CvMat*> intBeta;

		// Initialization
		// Reset memories
		memset(selected, 0, sizeof(bool) * N);
		// Computing the integral image of the original image
		vector<CvMat*> integralF;
		for (int i = 0; i < num_foregrounds; ++ i)
			integralF.push_back(calcIntegral(foreground[i]));
		vector<CvMat*> integralG;
		for (int i = 0; i < num_backgrounds; ++ i)
			integralG.push_back(calcIntegral(background[i]));

		/// Representative features from the dictionary
		int sampN = indexOfCentBases.size();

		int * progress = new int[N];

#define __VER2
#ifdef __VER2
		// Process of initializing parameters and get the first choice
		int tempL = 0;
		double maxFunc = -inf;
		double curFunc;
		for(int i = 0; i < N; ++ i) {
			progress[i] = 0;
			for (int j = 0; j < num_foregrounds; ++ j)
				BF[j][i] = listOfBases[i].innerProd(integralF[j]);
			for (int j = 0; j < num_backgrounds; ++ j)
				BG[j][i] = listOfBases[i].innerProd(integralG[j]);
			D[i] = 1;

			curFunc = 0;
			for (int j = 0; j < num_foregrounds; ++ j)
				curFunc += BF[j][i] * BF[j][i] / num_foregrounds;
			for (int j = 0; j < num_backgrounds; ++ j)
				curFunc -= BG[j][i] * BG[j][i] / num_backgrounds * lambda;

			if (curFunc > maxFunc) {
				tempL = i;
				maxFunc = curFunc;
			}
		}
		L.push_back(tempL);
		selected[tempL] = true;
#else
		// Process of initializing parameters and get the first choice
		int tempL = 0;
		double maxFunc = -inf;
		double curFunc;
		for(int i = 0; i < sampN; ++ i) {
			progress[i] = 0;
			for (int j = 0; j < num_foregrounds; ++ j)
				BF[j][i] = listOfBases[i].innerProd(integralF[j]);
			for (int j = 0; j < num_backgrounds; ++ j)
				BG[j][i] = listOfBases[i].innerProd(integralG[j]);
			D[i] = 1;

			curFunc = 0;
			for (int j = 0; j < num_foregrounds; ++ j)
				curFunc += BF[j][i] * BF[j][i] / num_foregrounds;
			for (int j = 0; j < num_backgrounds; ++ j)
				curFunc -= BG[j][i] * BG[j][i] / num_backgrounds * lambda;

			if (curFunc > maxFunc) {
				tempL = i;
				maxFunc = curFunc;
			}
		}
		L.push_back(tempL);
		selected[tempL] = true;
#endif

		// Computing parameters
		Phi.push_back(listOfBases[L[0]].toMatrix(height, width));
		intPhi.push_back(calcIntegral(Phi[0]));
		Beta.push_back(cvCloneMat(Phi[0]));
		intBeta.push_back(calcIntegral(Beta[0]));
		for (int i = 0; i < num_foregrounds; ++ i) {
			Coeff_fg.push_back(vector<double>());
			Coeff_fg[i].push_back(BF[i][L[0]]);
		}
		for (int i = 0; i < num_backgrounds; ++ i) {
			Coeff_bg.push_back(vector<double>());
			Coeff_bg[i].push_back(BG[i][L[0]]);
		}

		//Coeff.push_back(BF[L[0]]);
		for (int i = 0; i < num_foregrounds; ++ i)
			innerPhiF[i].push_back(innerProd(Phi[0], foreground[i]));
		for (int i = 0; i < num_backgrounds; ++ i)
			innerPhiG[i].push_back(innerProd(Phi[0], background[i]));

		int curtime = clock();

		double * sampScore = new double[sampN];
		double * bounds = new double[sampN];

		double vardelta = sqrt(2. - 2. * inprod_limit);

		K = 1;
		// Total Time Complexity is O(KN + K^2WH) = O(KWH(WH+K))
		while (K < MaxNumBases) {
			/**
			xmlDocPtr doc = xmlNewDoc(BAD_CAST"1.0");
			xmlNodePtr root_node = xmlNewNode(NULL, BAD_CAST"staterror");
			xmlDocSetRootElement(doc, root_node);
			static char buf[255];*/

			// Selecting the best fit Base
			maxFunc = -inf;
			tempL = -1;
			double tempPhiBase;
			int optIndex;

			for(int i = 0; i < sampN; ++ i) {
				//if (!selected[i]) {
				int ii = indexOfCentBases[i];
				progress[ii] = K;
				tempPhiBase = listOfBases[ii].innerProd(intPhi[K - 1]);
				for (int j = 0; j < num_foregrounds; ++ j)
					BF[j][ii] -= innerPhiF[j][K - 1] * tempPhiBase;
				for (int j = 0; j < num_backgrounds; ++ j)
					BG[j][ii] -= innerPhiG[j][K - 1] * tempPhiBase;
				D[ii] -= tempPhiBase * tempPhiBase;

				curFunc = 0;
				for (int j = 0; j < num_foregrounds; ++ j)
					curFunc += BF[j][ii] * BF[j][ii] / num_foregrounds;
				for (int j = 0; j < num_backgrounds; ++ j)
					curFunc -= BG[j][ii] * BG[j][ii] / num_backgrounds * lambda;

				curFunc = fabs(curFunc) < sqeps ? 0 : (curFunc / D[ii]);

				sampScore[i] = curFunc;

				/**
				xmlNodePtr node = xmlNewNode(NULL, BAD_CAST "base");
				xmlAddChild(root_node, node);
				xmlNewProp(node, BAD_CAST "index", BAD_CAST itoa(i, buf, 255));
				xmlNewProp(node, BAD_CAST "func", BAD_CAST double2string(curFunc, buf, 255)); */

				if (curFunc > maxFunc) {
					maxFunc = curFunc;
					tempL = ii;
					optIndex = i;
				}
			}

			/**
			/// Save XML documents
			//Dumping document to stdio or file
			xmlSaveFormatFileEnc((string("exp/") + itoa(K, buf, 255) + string(".xml")).c_str(), doc, "UTF-8", 1);
			xmlFreeDoc(doc);
			xmlCleanupParser();
			xmlMemoryDump();//debug memory for regression tests */

			if (tempL == -1) break;

			//double templ_thresh = maxFunc - fabs(maxFunc) * ratio;

			for (int ind = 0; ind < sampN; ++ ind)
				if (sampScore[ind] > maxFunc - fabs(maxFunc) * ratio) {

					for (int i = 0; i < (int)memberOfGroup[ind].size(); ++ i) {
						int ii = memberOfGroup[ind][i];

						for (int k = progress[ii]; k < K; ++ k) {
							tempPhiBase = listOfBases[ii].innerProd(intPhi[k]);
							for (int j = 0; j < num_foregrounds; ++ j)
								BF[j][ii] -= innerPhiF[j][k] * tempPhiBase;
							for (int j = 0; j < num_backgrounds; ++ j)
								BG[j][ii] -= innerPhiG[j][k] * tempPhiBase;
							D[ii] -= tempPhiBase * tempPhiBase;
						}

						progress[ii] = K;

						curFunc = 0;
						for (int j = 0; j < num_foregrounds; ++ j)
							curFunc += BF[j][ii] * BF[j][ii] / num_foregrounds;
						for (int j = 0; j < num_backgrounds; ++ j)
							curFunc -= BG[j][ii] * BG[j][ii] / num_backgrounds * lambda;

						curFunc = fabs(curFunc) < sqeps ? 0 : (curFunc / D[ii]);

						if (curFunc > maxFunc) {
							maxFunc = curFunc;
							tempL = ii;
						}

					}
				}

				L.push_back(tempL);
				selected[tempL] = true;

				// Computing parameters
				CvMat * tempMat = calcPhi(L[K], Phi, intPhi);       // O(KWH)
				Phi.push_back(calcMatDiv(tempMat, sqrt(D[L[K]])));  // O(WH)
				Beta.push_back(calcMatDiv(tempMat, D[L[K]]));       // O(WH)
				intBeta.push_back(calcIntegral(Beta[K]));           // O(WH)
				cvReleaseMat(&tempMat);
				intPhi.push_back(calcIntegral(Phi[K]));             // O(WH)
				for (int i = 0; i < num_foregrounds; ++ i)
					innerPhiF[i].push_back(innerProd(Phi[K], foreground[i]));
				for (int i = 0; i < num_backgrounds; ++ i)
					innerPhiG[i].push_back(innerProd(Phi[K], background[i]));

				// Next iteration
				++ K;
		}

		lyon::settime((clock() - curtime) / (double)CLOCKS_PER_SEC);

		// Copy Answers
		this->indexOfBases.clear();
		for(int i = 0; i < (int)L.size(); ++ i)
			this->indexOfBases.push_back(listOfBases[L[i]]);
		CvMat * intTarget = calcIntegral(target);
		this->coefMatrix();
		this->coeffOfBases = compCoeff(intTarget);

		//demo();

		// Release Memories
		delete [] progress;
		delete [] D;
		delete [] selected;
		for (int i = 0; i < num_foregrounds; ++ i) {
			delete [] BF[i];
		}
		for (int i = 0; i < num_backgrounds; ++ i) {
			delete [] BG[i];
		}
		for(int i = 0; i < (int)Phi.size(); ++ i) {
			cvReleaseMat(&Phi[i]);
			cvReleaseMat(&intPhi[i]);
			cvReleaseMat(&Beta[i]);
			cvReleaseMat(&intBeta[i]);
		}
}




void DnbsClassGamma::computeDNBS(int MaxNumBases,
								 int MaxNumForegrounds,
								 int MaxNumBackgrounds,
								 IplImage * frame,
								 CvRect rect,
								 int margin_threshold_x,
								 int margin_threshold_y,
								 int distance_threshold,
								 double lambda,
								 int method,
								 int samp_width,
								 int samp_height,
								 double prune_ratio,
								 CvMat ** previous_foregrounds,
								 double coherence) 
{
	// Ang Li @ Nanjing Univ. Jan.27, 2010.

	// First compute the original NBS representation
	//importTarget(target);
	computeNBS_v1_cluster_v2(30, 0.5);
	//NaiveNbsClass::computeNBS(30);
	//printf("Original Reconstruction error: %.6lf\n", reconerr());

	/*
	if (MaxNumForegrounds == 1 && MaxNumBackgrounds == 0) {
	indexOfBases = naivenbs.indexOfBases;
	coefMatrix();
	return ;
	}
	*/

	// Second update the NBS by sampling background templates

	// Get the window of neighborhood
	int margin_x = margin_threshold_x;//max(margin_threshold_x, rect.width);
	int margin_y = margin_threshold_y;//max(margin_threshold_y, rect.height);

	// Construct the SSD distribution of the neighborhood
	double fnorm = l2norm();
	//double denom_beta = 256. * 256. * rect.width * rect.height;
	CvMat * intframe = lyon::calcIntegral(frame, 0, frame->width, 0, frame->height);
	CvMat * int2frame = lyon::calcSqrIntegral(frame, 0, frame->width, 0, frame->height);

	int x_min = max(rect.x - margin_x, 0);
	int x_max = min(rect.x + margin_x, frame->width - imagW);
	int y_min = max(rect.y - margin_y, 0);
	int y_max = min(rect.y + margin_y, frame->height - imagH);
	CvSize neighborSize = cvSize(x_max - x_min, y_max - y_min);
	double ** weight = (double**)malloc(sizeof(double*) * neighborSize.width);
	for (int i = 0; i < neighborSize.width; ++ i) {
		weight[i] = (double*)malloc(sizeof(double) * neighborSize.height);
		memset(weight[i], 0, sizeof(double) * neighborSize.height);
	}

	vector< pair<double, int> > arr;
	vector< CvRect > pnt;

	// Calculate SSD between neighbor patches and the target
	// approximated by previous constructed NBS
	// weight[x - x_min][y - y_min] = SSD(x,y)
	for (int x = x_min; x < x_max; ++ x) {
		for (int y = y_min; y < y_max; ++ y) {
			double gnorm
				= lyon::calcIntSum(int2frame, x, x + imagW, y, y + imagH);
			double iproduct
				= innerProduct(intframe, x, x + imagW, y, y + imagH);
			double weight_of_patch = (gnorm - 2. * iproduct + fnorm);
			weight[x - x_min][y - y_min] = weight_of_patch;
			pnt.push_back(cvRect(x, y, target->width, target->height));
			arr.push_back(make_pair(weight_of_patch, pnt.size() - 1));
		}
	}

	// Sample background templates in light of the SSD distribution
	// Sort 'arr' in ascending order of weights
	sort(arr.begin(), arr.end());

	// Set sample vectors and their weights
	vector < CvRect > samples;

#define _BG_SAMPLE_ 1
#if (_BG_SAMPLE_ == 1)
 	sample_backgrounds(
 		rect, MaxNumBackgrounds, distance_threshold,
 		x_min, x_max, y_min, y_max, pnt, arr, weight, samples);
#elif (_BG_SAMPLE_ == 2)
	sample_backgrounds_supress(
		rect, MaxNumBackgrounds, distance_threshold,
		x_min, x_max, y_min, y_max, pnt, arr, weight, samples, 0.7);
#endif

	printf ("Number of sampled backgrounds = %d\n", samples.size());

	// bgsample records the sampled backgrounds
	bgsample.clear();
	for (int i = 0; i < (int)samples.size(); ++ i)
		bgsample.push_back(samples[i]);

	// Collect positive and negative samples
	vector<CvMat*> foreground;
	vector<CvMat*> background;
	//foreground.push_back(target);
	for (int i = 0; i < MaxNumForegrounds /*- 1*/; ++ i)
		if (previous_foregrounds[i] != NULL)
			foreground.push_back(previous_foregrounds[i]);
	for (int i = 0; i < (int)samples.size(); ++ i) {
		IplImage * patch
			= lyon::selectedPatch(frame, samples[i]);
		background.push_back(image2mat(patch));
		cvReleaseImage(&patch);
	}

	// Re-compute the NBS of foreground that discriminates the backgrounds
	switch (method) {
		case 0:
			computeDNBS_v1(MaxNumBases, lambda, foreground, background);
			break;
		case 1:
			computeDNBS_v2(MaxNumBases, frame, lambda, samples);
			break;
		case 2:
			computeDNBS_v1_pruning(
				MaxNumBases, lambda, foreground, background, samp_width, samp_height);
			break;
		case 3:
			computeDNBS_v1_downsample(
				MaxNumBases, lambda, foreground, background, samp_width, samp_height);
			break;
		case 4:
			computeDNBS_v1_pruning_2nd(
				MaxNumBases, lambda, foreground, background, samp_width, samp_height);
			break;
		case 5:
			computeDNBS_v1_branchbound(
				MaxNumBases, lambda, foreground, background, inprod_limit);
			break;
		case 6:
			computeDNBS_v1_cluster(
				MaxNumBases, lambda, foreground, background, prune_ratio);
			break;
		case 7:
			computeDNBS_v1_branchbound_latest(
				MaxNumBases, lambda, foreground, background);
			break;
		case 8: { // downsample and then work on the generated new feature set
			vector<int> featIndex = genDownScaleFeatures(samp_width, samp_height);
			featIndex = computeDNBS_v1_customFeatSet(
				MaxNumBases, lambda, foreground, background, featIndex);
			featIndex = genDiffuseFeatures(featIndex, 0.75);
			computeDNBS_v1_customFeatSet(
				MaxNumBases, lambda, foreground, background, featIndex);
			break;
				}
		case 9: // Using clusters, re-score the first rank cluster
			computeDNBS_v1_cluster_v2(
				MaxNumBases, lambda, foreground, background, prune_ratio);
			break;
		case 10:
			computeDNBS_v2_cluster(
				MaxNumBases, lambda, foreground, background, prune_ratio);
			break;
		default:
			break;
	}

	// Back-up
	for (int i = 0; i < (int)fgrdSample.size(); ++ i)
		cvReleaseMat(&fgrdSample[i]);
	for (int i = 0; i < (int)bgrdSample.size(); ++ i)
		cvReleaseMat(&bgrdSample[i]);
	fgrdSample.clear();
	bgrdSample.clear();
	for (int i = 0; i < (int)foreground.size(); ++ i)
		fgrdSample.push_back(cvCloneMat(foreground[i]));
	for (int i = 0; i < (int)background.size(); ++ i)
		bgrdSample.push_back(cvCloneMat(background[i]));

	// Release arrays and images
	for (int i = 0; i < (int)background.size(); ++ i)
		cvReleaseMat(&background[i]);
	//printf("Neighbor Size = (%d,%d)\n", neighborSize.width, neighborSize.height);
	for (int i = 0; i < neighborSize.width; ++ i)
		free(weight[i]);
	free(weight);
	cvReleaseMat(&intframe);
	cvReleaseMat(&int2frame);

}

bool DnbsClassGamma::genHaarBases_innProdFast(double minInnProd) {
	// This method implements the fast feature clustering in accord
	// with the given limitation of mutual inner products, i.e. the
	// inner product between each pair of features within the same
	// cluster should not below 'minInnProd'.
	// The data members are clarified:
	// 1. 'double inprod_limit' stores the limitation of products
	// 2. 'vector<int> indexOfCentBases' stores the index of the center
	//    basis in each of the clusters.
	// 3. 'vector< vector<int> > memberOfGroup' stores the index of
	//    features in each of the clusters.
	// This version of feature grouping is faster, because optimized
	// by considering the connection between common areas and inner
	// products for the Haar-like rectangular feature set.
	// Ang Li @ Nanjing University. March 2011.

	// Initialization: clear old data
	indexOfCentBases.clear();
	for (int i = 0; i < (int)memberOfGroup.size(); ++ i)
		memberOfGroup[i].clear();
	memberOfGroup.clear();

	// Initialization: assign new value to inprod_limit
	inprod_limit = minInnProd;

	// Efficient checking according to common areas.

	int N = listOfBases.size();
	int remainN = N;
	int * locIndex = new int[N];
	int * remainIndex = new int[N];
	for (int i = 0; i < N; ++ i)
		locIndex[i] = remainIndex[i] = i;

	// 	bool * visited = new bool[N];
	// 	memset(visited, false, sizeof(bool) * N);

	int cnttime = clock();
	//srand(time(NULL));

	while (remainN > 0) {
		//printf("remainN = %d\n", remainN);
		//printf("# cluster = %d, # remain = %d\n", indexOfCentBases.size(), remainN);

		int ind = (int)(rand() / (double)RAND_MAX * remainN);
		//assert(ind < remainN);
#ifdef __DEBUG
		if (indexOfCentBases.size() < 10) {
			printf("remainN = %d\n", remainN);
			printf("ind = %d\n", ind);
		}
#endif
		//visited[remainIndex[ind]] = true;
		int cind = remainIndex[ind];
		indexOfCentBases.push_back(cind);
		memberOfGroup.push_back(vector<int>());
		vector<int> & locMemOfGroup = memberOfGroup[memberOfGroup.size() - 1];
		// clear up nearby features
		swap(remainIndex[ind], remainIndex[remainN - 1]);
		swap(locIndex[remainIndex[ind]], locIndex[remainIndex[remainN - 1]]);
		-- remainN;
		HaarBase & locHaarBase = listOfBases[cind];
		int oW = locHaarBase.w;
		int oH = locHaarBase.h;

		double tmplthresh = oW * oH * minInnProd * minInnProd;
		//system("pause");

		int leftlimit, rightlimit, toplimit, botlimit;

		for (int lx = 0; lx < oW; ++ lx)
			for (int rx = lx; rx < oW; ++ rx) {
				int pw = rx - lx + 1;
				int miniH = (int)ceil(tmplthresh / pw - eps);
				for (int uy = 0; uy < oH; ++ uy)
					for (int by = uy + miniH - 1; by < oH; ++ by) {
						int ph = by - uy + 1;
						int totlimit = (int)((double)ph * ph * pw * pw / tmplthresh);
						leftlimit = (lx == 0) ? min(locHaarBase.wtl, totlimit / ph - pw) : 0;
						for (int extleft = 0; extleft <= leftlimit; ++ extleft) {
							if (rx == oW - 1) {
								rightlimit = min(imagW - locHaarBase.wtl - locHaarBase.w, totlimit / ph - pw - extleft);
							} else rightlimit = 0;
							for (int extright = 0; extright <= rightlimit; ++ extright) {
								if (uy == 0) {
									toplimit = totlimit / (pw + extleft + extright) - ph;
									toplimit = min(toplimit, locHaarBase.htl);
								} else toplimit = 0;
								for (int exttop = 0; exttop <= toplimit; ++ exttop) {
									if (by == oH - 1) {
										botlimit = (int)(totlimit / (pw + extleft + extright) - ph - exttop);
										botlimit = min(botlimit, imagH - locHaarBase.htl - locHaarBase.h);
									} else botlimit = 0;
									for (int extbot = 0; extbot <= botlimit; ++ extbot) {
										// this basis should be included in the current cluster.
										// (extleft, extright, exttop, extbot)
										//printf("limit: (%d, %d, %d, %d)\n", leftlimit, rightlimit, toplimit, botlimit);

										//printf("(%d, %d, %d, %d)\n", extleft, extright, exttop, extbot);
										int ii = calc_index_haarbase(
											locHaarBase.htl + uy - exttop,
											locHaarBase.wtl + lx - extleft,
											ph + exttop + extbot,
											pw + extleft + extright,
											imagH, imagW);
										//if (visited[ii]) continue;
										if (locIndex[ii] >= remainN) continue;
										// 										if (locHaarBase.innerProd(listOfBases[ii]) < minInnProd - eps) {
										// 											printf("Std Inner Prod = %.3lf\n", locHaarBase.innerProd(listOfBases[ii]));
										// 											printf("Quick Inner Prod = %.3lf\n", ph * pw / sqrt((double)oW * oH * (ph + exttop + extbot) * (pw + extleft + extright)));
										// 											system("pause");
										// 											continue;
										// 										}
										//visited[ii] = true;
										locMemOfGroup.push_back(ii);
										swap(locIndex[ii], locIndex[remainIndex[remainN - 1]]);
										swap(remainIndex[locIndex[ii]], remainIndex[locIndex[remainIndex[remainN - 1]]]);
										-- remainN;
									}
								}
							}
						}
					}
			}

	}

	printf("Time consumed: %.3lf\n", (clock() - cnttime) / (double)CLOCKS_PER_SEC);
	printf("Hierarchy generated: %d clusters / %d.\n", indexOfCentBases.size(), N);

	delete [] locIndex;
	delete [] remainIndex;

	return true;
}

void DnbsClassGamma::computeDNBS_v1_cluster_v2(
	int MaxNumBases,
	double lambda,
	vector<CvMat*> foreground,
	vector<CvMat*> background,
	double ratio) {

		// Supposing frame is gray-level image with 3 channels

		// Supposing 'GenHaarBases_innerprod' completed.

		// Record the width and height
		int width = foreground[0]->width;
		int height = foreground[0]->height;

		int num_foregrounds = foreground.size();
		int num_backgrounds = background.size();

		// Definitions
		// K is the number of selected bases
		int K;
		// N is the total number of candidate bases
		int N = listOfBases.size();

		// B and D are two arrays of parameters in deciding the
		// optimized index of the next to-be-selected binary base.
		// i.e. selecting the one that maximising E_n(=|BF_n||BG_n|/D_n)
		vector<NBS_MAT_TYPE*> BF;
		for (int i = 0; i < num_foregrounds; ++ i) {
			NBS_MAT_TYPE * tempBF
				= new NBS_MAT_TYPE[N];
			BF.push_back(tempBF);
		}

		vector<NBS_MAT_TYPE*> BG;
		for (int i = 0; i < num_backgrounds; ++ i) {
			NBS_MAT_TYPE * tempBG
				= new NBS_MAT_TYPE[N];
			BG.push_back(tempBG);
		}
		NBS_MAT_TYPE * D = new NBS_MAT_TYPE[N];

		// Indicating whether the base is selected to form the subspace
		bool * selected = new bool[N];

		// Index of each selected binary bases
		vector<int> L;
		// Projection Coefficients
		vector<NBS_MAT_TYPE> Coeff;
		vector< vector<double> > Coeff_fg, Coeff_bg;

		// Phi's and their integral image
		vector<CvMat*> Phi;
		vector<CvMat*> intPhi;

		// Inner product of each Phi[K] and image F(i.e. matPatch)
		vector< vector<NBS_MAT_TYPE> > innerPhiF;
		for (int i = 0; i < num_foregrounds; ++ i) 
			innerPhiF.push_back(vector<NBS_MAT_TYPE>());
		// Inner product of each Phi[K] and image G
		vector< vector<NBS_MAT_TYPE> > innerPhiG;
		for (int i = 0; i < num_backgrounds; ++ i)
			innerPhiG.push_back(vector<NBS_MAT_TYPE>());

		// Beta's and their integral image
		vector<CvMat*> Beta;
		vector<CvMat*> intBeta;

		// Initialization
		// Reset memories
		memset(selected, 0, sizeof(bool) * N);
		// Computing the integral image of the original image
		vector<CvMat*> integralF;
		for (int i = 0; i < num_foregrounds; ++ i)
			integralF.push_back(calcIntegral(foreground[i]));
		vector<CvMat*> integralG;
		for (int i = 0; i < num_backgrounds; ++ i)
			integralG.push_back(calcIntegral(background[i]));

		/// Representative features from the dictionary
		int sampN = indexOfCentBases.size();

		int * progress = new int[N];

		// Process of initializing parameters and get the first choice
		int tempL = 0;
		double maxFunc = -inf;
		double curFunc;
		for(int i = 0; i < N; ++ i) {
			progress[i] = 0;
			for (int j = 0; j < num_foregrounds; ++ j)
				BF[j][i] = listOfBases[i].innerProd(integralF[j]);
			for (int j = 0; j < num_backgrounds; ++ j)
				BG[j][i] = listOfBases[i].innerProd(integralG[j]);
			D[i] = 1;

			curFunc = 0;
			for (int j = 0; j < num_foregrounds; ++ j)
				curFunc += BF[j][i] * BF[j][i] / num_foregrounds;
			for (int j = 0; j < num_backgrounds; ++ j)
				curFunc -= BG[j][i] * BG[j][i] / num_backgrounds * lambda;

			if (curFunc > maxFunc) {
				tempL = i;
				maxFunc = curFunc;
			}
		}
		L.push_back(tempL);
		selected[tempL] = true;

		// Computing parameters
		Phi.push_back(listOfBases[L[0]].toMatrix(height, width));
		intPhi.push_back(calcIntegral(Phi[0]));
		Beta.push_back(cvCloneMat(Phi[0]));
		intBeta.push_back(calcIntegral(Beta[0]));
		for (int i = 0; i < num_foregrounds; ++ i) {
			Coeff_fg.push_back(vector<double>());
			Coeff_fg[i].push_back(BF[i][L[0]]);
		}
		for (int i = 0; i < num_backgrounds; ++ i) {
			Coeff_bg.push_back(vector<double>());
			Coeff_bg[i].push_back(BG[i][L[0]]);
		}

		//Coeff.push_back(BF[L[0]]);
		for (int i = 0; i < num_foregrounds; ++ i)
			innerPhiF[i].push_back(innerProd(Phi[0], foreground[i]));
		for (int i = 0; i < num_backgrounds; ++ i)
			innerPhiG[i].push_back(innerProd(Phi[0], background[i]));

		int curtime = clock();

		double * sampScore = new double[sampN];
		double * bounds = new double[sampN];

		double vardelta = sqrt(2. - 2. * inprod_limit);

		//#define _BOUND_PLOT
#ifdef _BOUND_PLOT
		char filename[225];
		sprintf(filename, "bound_stat_%.3lf.txt", inprod_limit);
		FILE * fout = fopen(filename, "w+");
#endif

		K = 1;
		// Total Time Complexity is O(KN + K^2WH) = O(KWH(WH+K))
		while (K < MaxNumBases) {
			/**
			xmlDocPtr doc = xmlNewDoc(BAD_CAST"1.0");
			xmlNodePtr root_node = xmlNewNode(NULL, BAD_CAST"staterror");
			xmlDocSetRootElement(doc, root_node);
			static char buf[255];*/

			// Selecting the best fit Base
			maxFunc = -inf;
			tempL = -1;
			double tempPhiBase;
			int optIndex;

			for(int i = 0; i < sampN; ++ i) {
				//if (!selected[i]) {
				int ii = indexOfCentBases[i];
				progress[ii] = K;
				tempPhiBase = listOfBases[ii].innerProd(intPhi[K - 1]);
				for (int j = 0; j < num_foregrounds; ++ j)
					BF[j][ii] -= innerPhiF[j][K - 1] * tempPhiBase;
				for (int j = 0; j < num_backgrounds; ++ j)
					BG[j][ii] -= innerPhiG[j][K - 1] * tempPhiBase;
				D[ii] -= tempPhiBase * tempPhiBase;

				curFunc = 0;
				for (int j = 0; j < num_foregrounds; ++ j)
					curFunc += BF[j][ii] * BF[j][ii] / num_foregrounds;
				for (int j = 0; j < num_backgrounds; ++ j)
					curFunc -= BG[j][ii] * BG[j][ii] / num_backgrounds * lambda;

				curFunc = fabs(D[ii]) < sqeps ? 0 : (curFunc / D[ii]);

				sampScore[i] = curFunc;

				/**
				xmlNodePtr node = xmlNewNode(NULL, BAD_CAST "base");
				xmlAddChild(root_node, node);
				xmlNewProp(node, BAD_CAST "index", BAD_CAST itoa(i, buf, 255));
				xmlNewProp(node, BAD_CAST "func", BAD_CAST double2string(curFunc, buf, 255)); */

				if (curFunc > maxFunc) {
					maxFunc = curFunc;
					tempL = ii;
					optIndex = i;
				}
			}

			/**
			/// Save XML documents
			//Dumping document to stdio or file
			xmlSaveFormatFileEnc((string("exp/") + itoa(K, buf, 255) + string(".xml")).c_str(), doc, "UTF-8", 1);
			xmlFreeDoc(doc);
			xmlCleanupParser();
			xmlMemoryDump();//debug memory for regression tests */

			if (tempL == -1) break;

			//double templ_thresh = maxFunc - fabs(maxFunc) * ratio;

			for (int ind = 0; ind < sampN; ++ ind)
#ifdef _BOUND_PLOT
			{
#else
				if (sampScore[ind] > maxFunc - fabs(maxFunc) * ratio) {
#endif

					HaarBase & locHaarBase = listOfBases[indexOfCentBases[ind]];
					int oW = locHaarBase.w;
					int oH = locHaarBase.h;

					double tmplthresh = oW * oH * inprod_limit * inprod_limit;

					int leftlimit, rightlimit, toplimit, botlimit;

#ifdef _BOUND_PLOT
					double maxbound = 0;
#endif

					for (int lx = 0; lx < oW; ++ lx)
						for (int rx = lx; rx < oW; ++ rx) {
							int pw = rx - lx + 1;
							int miniH = (int)ceil(tmplthresh / pw - eps);
							for (int uy = 0; uy < oH; ++ uy)
								for (int by = uy + miniH - 1; by < oH; ++ by) {
									int ph = by - uy + 1;
									int totlimit = (int)((double)ph * ph * pw * pw / tmplthresh);
									leftlimit = (lx == 0) ? min(locHaarBase.wtl, totlimit / ph - pw) : 0;
									for (int extleft = 0; extleft <= leftlimit; ++ extleft) {
										if (rx == oW - 1) {
											rightlimit = min(imagW - locHaarBase.wtl - locHaarBase.w, totlimit / ph - pw - extleft);
										} else rightlimit = 0;
										for (int extright = 0; extright <= rightlimit; ++ extright) {
											if (uy == 0) {
												toplimit = totlimit / (pw + extleft + extright) - ph;
												toplimit = min(toplimit, locHaarBase.htl);
											} else toplimit = 0;
											for (int exttop = 0; exttop <= toplimit; ++ exttop) {
												if (by == oH - 1) {
													botlimit = (int)(totlimit / (pw + extleft + extright) - ph - exttop);
													botlimit = min(botlimit, imagH - locHaarBase.htl - locHaarBase.h);
												} else botlimit = 0;
												for (int extbot = 0; extbot <= botlimit; ++ extbot) {
													// this basis should be included in the current cluster.
													// (extleft, extright, exttop, extbot)
													int ii = calc_index_haarbase(
														locHaarBase.htl + uy - exttop,
														locHaarBase.wtl + lx - extleft,
														ph + exttop + extbot,
														pw + extleft + extright,
														imagH, imagW);

#ifdef _BOUND_PLOT

#else
													if (progress[ii] == K) continue;
#endif

													// check basis-ii to update the optimal value
													for (int k = progress[ii]; k < K; ++ k) {
														tempPhiBase = listOfBases[ii].innerProd(intPhi[k]);
														for (int j = 0; j < num_foregrounds; ++ j)
															BF[j][ii] -= innerPhiF[j][k] * tempPhiBase;
														for (int j = 0; j < num_backgrounds; ++ j)
															BG[j][ii] -= innerPhiG[j][k] * tempPhiBase;
														D[ii] -= tempPhiBase * tempPhiBase;
													}

													progress[ii] = K;

													curFunc = 0;
													for (int j = 0; j < num_foregrounds; ++ j)
														curFunc += BF[j][ii] * BF[j][ii] / num_foregrounds;
													for (int j = 0; j < num_backgrounds; ++ j)
														curFunc -= BG[j][ii] * BG[j][ii] / num_backgrounds * lambda;

													curFunc = fabs(D[ii]) < sqeps ? 0 : (curFunc / D[ii]);

#ifdef _BOUND_PLOT
													if (curFunc - sampScore[ind] > maxbound) {
														maxbound = curFunc - sampScore[ind];
													}
#endif

													if (curFunc > maxFunc) {
														maxFunc = curFunc;
														tempL = ii;
													}
												}
											}
										}
									}
								}
						}

#ifdef _BOUND_PLOT
						// K, Bound, SampScore, D
						fprintf(fout, "%d\t%lf\t%lf\t%lf\n", K, maxbound, sampScore[ind], D[indexOfCentBases[ind]]);
#endif

				}

				// 				// check
				// 				int cnt = 0;
				// 				for (int i = 0; i < N; ++ i)
				// 					if (progress[i] != K) {
				// 						++ cnt;
				// 					}
				// 					printf("#unvisited / N = %d / %d\n", cnt, N);

				L.push_back(tempL);
				selected[tempL] = true;

				// Computing parameters
				CvMat * tempMat = calcPhi(L[K], Phi, intPhi);       // O(KWH)
				Phi.push_back(calcMatDiv(tempMat, sqrt(D[L[K]])));  // O(WH)
				Beta.push_back(calcMatDiv(tempMat, D[L[K]]));       // O(WH)
				intBeta.push_back(calcIntegral(Beta[K]));           // O(WH)
				cvReleaseMat(&tempMat);
				intPhi.push_back(calcIntegral(Phi[K]));             // O(WH)
				for (int i = 0; i < num_foregrounds; ++ i)
					innerPhiF[i].push_back(innerProd(Phi[K], foreground[i]));
				for (int i = 0; i < num_backgrounds; ++ i)
					innerPhiG[i].push_back(innerProd(Phi[K], background[i]));

				// Next iteration
				++ K;
		}

#ifdef _BOUND_PLOT
		fclose(fout);
#endif

		lyon::settime((clock() - curtime) / (double)CLOCKS_PER_SEC);

		// Copy Answers
		this->indexOfBases.clear();
		for(int i = 0; i < (int)L.size(); ++ i)
			this->indexOfBases.push_back(listOfBases[L[i]]);
		CvMat * intTarget = calcIntegral(target);
		this->coefMatrix();
		this->coeffOfBases = compCoeff(intTarget);

		//demo();

		// Release Memories
		delete [] progress;
		delete [] D;
		delete [] selected;
		for (int i = 0; i < num_foregrounds; ++ i) {
			delete [] BF[i];
		}
		for (int i = 0; i < num_backgrounds; ++ i) {
			delete [] BG[i];
		}
		for(int i = 0; i < (int)Phi.size(); ++ i) {
			cvReleaseMat(&Phi[i]);
			cvReleaseMat(&intPhi[i]);
			cvReleaseMat(&Beta[i]);
			cvReleaseMat(&intBeta[i]);
		}
}

void DnbsClassGamma::computeDNBS_v2_cluster(
	int MaxNumBases,
	double lambda,
	vector<CvMat*> foreground,
	vector<CvMat*> background,
	double ratio) {

		// Supposing frame is gray-level image with 3 channels

		// Supposing 'GenHaarBases_innerprod' completed.

		// Record the width and height
		int width = foreground[0]->width;
		int height = foreground[0]->height;

		int num_foregrounds = foreground.size();
		int num_backgrounds = background.size();

		// Definitions
		// K is the number of selected bases
		int K;
		// N is the total number of candidate bases
		int N = listOfBases.size();

		// B and D are two arrays of parameters in deciding the
		// optimized index of the next to-be-selected binary base.
		// i.e. selecting the one that maximising E_n(=|BF_n||BG_n|/D_n)
		vector<NBS_MAT_TYPE*> BF;
		for (int i = 0; i < num_foregrounds; ++ i) {
			NBS_MAT_TYPE * tempBF
				= new NBS_MAT_TYPE[N];
			BF.push_back(tempBF);
		}

		vector<NBS_MAT_TYPE*> BG;
		for (int i = 0; i < num_backgrounds; ++ i) {
			NBS_MAT_TYPE * tempBG
				= new NBS_MAT_TYPE[N];
			BG.push_back(tempBG);
		}
		NBS_MAT_TYPE * D = new NBS_MAT_TYPE[N];

		// Indicating whether the base is selected to form the subspace
		bool * selected = new bool[N];

		// Index of each selected binary bases
		vector<int> L;
		// Projection Coefficients
		vector<NBS_MAT_TYPE> Coeff;
		vector< vector<double> > Coeff_fg, Coeff_bg;

		// Phi's and their integral image
		vector<CvMat*> Phi;
		vector<CvMat*> intPhi;

		// Inner product of each Phi[K] and image F(i.e. matPatch)
		vector< vector<NBS_MAT_TYPE> > innerPhiF;
		for (int i = 0; i < num_foregrounds; ++ i) 
			innerPhiF.push_back(vector<NBS_MAT_TYPE>());
		// Inner product of each Phi[K] and image G
		vector< vector<NBS_MAT_TYPE> > innerPhiB;
		for (int i = 0; i < num_backgrounds; ++ i)
			innerPhiB.push_back(vector<NBS_MAT_TYPE>());


		// Initialization
		// Reset memories
		memset(selected, 0, sizeof(bool) * N);
		// Computing the integral image of the original image
		vector<CvMat*> integralF;
		for (int i = 0; i < num_foregrounds; ++ i)
			integralF.push_back(calcIntegral(foreground[i]));
		vector<CvMat*> integralB;
		for (int i = 0; i < num_backgrounds; ++ i)
			integralB.push_back(calcIntegral(background[i]));

		// Computing the integral image of the original image
		vector<CvMat*> epsilon_f;
		for (int i = 0; i < num_foregrounds; ++ i)
			epsilon_f.push_back(cvCloneMat(foreground[i]));
		vector<CvMat*> epsilon_b;
		for (int i = 0; i < num_backgrounds; ++ i)
			epsilon_b.push_back(cvCloneMat(background[i]));
		vector<CvMat*> int_epsilon_f;
		for (int i = 0; i < num_foregrounds; ++ i)
			int_epsilon_f.push_back(cvCloneMat(integralF[i]));
		vector<CvMat*> int_epsilon_b;
		for (int i = 0; i < num_backgrounds; ++ i)
			int_epsilon_b.push_back(cvCloneMat(integralB[i]));

		/// Representative features from the dictionary
		int sampN = indexOfCentBases.size();
		double *sampScore = new double[sampN];

		double *Lvalue = new double[N];

		int * progress = new int[N];
		memset(progress, -1, sizeof(int) * N);

		// Selecting the best fit Base
        double maxFunc = -inf;
        int tempL = -1;

		for(int index = 0; index < N; ++ index) {
			//if (!selected[i]) {
			progress[index] = 0;
			double sum_pos = 0, sum_neg = 0;
			for (int j = 0; j < num_foregrounds; ++ j)
				sum_pos += lyon::sqr(listOfBases[index].innerProd(integralF[j]));
			for (int j = 0; j < num_backgrounds; ++ j)
				sum_neg += lyon::sqr(listOfBases[index].innerProd(integralB[j]));
			Lvalue[index] = sum_pos / num_foregrounds - (num_backgrounds ? lambda * sum_neg / num_backgrounds : 0);
			D[index] = 1;
			if (Lvalue[index] > maxFunc) {
				tempL = index;
				maxFunc = Lvalue[index];
			}
		}

#if 0
		for(int i = 0; i < sampN; ++ i) {
			//if (!selected[i]) {
			int index = indexOfCentBases[i];
			progress[index] = 0;
			double sum_pos = 0, sum_neg = 0;
			for (int j = 0; j < num_foregrounds; ++ j)
				sum_pos += lyon::sqr(listOfBases[index].innerProd(integralF[j]));
			for (int j = 0; j < num_backgrounds; ++ j)
				sum_neg += lyon::sqr(listOfBases[index].innerProd(integralB[j]));
			Lvalue[index] = sum_pos / num_foregrounds - lambda * sum_neg / num_backgrounds;
			D[index] = 1;
			if (Lvalue[index] > maxFunc) {
				tempL = index;
				maxFunc = Lvalue[index];
			}
			sampScore[i] = Lvalue[index];
		}

		for (int ind = 0; ind < sampN; ++ ind)
			if (sampScore[ind] > maxFunc - fabs(maxFunc) * ratio) {

				HaarBase & locHaarBase = listOfBases[indexOfCentBases[ind]];
				int oW = locHaarBase.w;
				int oH = locHaarBase.h;

				double tmplthresh = oW * oH * inprod_limit * inprod_limit;

				int leftlimit, rightlimit, toplimit, botlimit;

				for (int lx = 0; lx < oW; ++ lx)
					for (int rx = lx; rx < oW; ++ rx) {
						int pw = rx - lx + 1;
						int miniH = (int)ceil(tmplthresh / pw - eps);
						for (int uy = 0; uy < oH; ++ uy)
							for (int by = uy + miniH - 1; by < oH; ++ by) {
								int ph = by - uy + 1;
								int totlimit = (int)((double)ph * ph * pw * pw / tmplthresh);
								leftlimit = (lx == 0) ? min(locHaarBase.wtl, totlimit / ph - pw) : 0;
								for (int extleft = 0; extleft <= leftlimit; ++ extleft) {
									if (rx == oW - 1) {
										rightlimit = min(imagW - locHaarBase.wtl - locHaarBase.w, totlimit / ph - pw - extleft);
									} else rightlimit = 0;
									for (int extright = 0; extright <= rightlimit; ++ extright) {
										if (uy == 0) {
											toplimit = totlimit / (pw + extleft + extright) - ph;
											toplimit = min(toplimit, locHaarBase.htl);
										} else toplimit = 0;
										for (int exttop = 0; exttop <= toplimit; ++ exttop) {
											if (by == oH - 1) {
												botlimit = (int)(totlimit / (pw + extleft + extright) - ph - exttop);
												botlimit = min(botlimit, imagH - locHaarBase.htl - locHaarBase.h);
											} else botlimit = 0;
											for (int extbot = 0; extbot <= botlimit; ++ extbot) {
												// this basis should be included in the current cluster.
												// (extleft, extright, exttop, extbot)
												int ii = calc_index_haarbase(
													locHaarBase.htl + uy - exttop,
													locHaarBase.wtl + lx - extleft,
													ph + exttop + extbot,
													pw + extleft + extright,
													imagH, imagW);

												if (progress[ii] == 0) continue;

												if (progress[ii] == -1) {
													HaarBase & cbase = listOfBases[ii];
													progress[ii] = 0;
													double sum_pos = 0, sum_neg = 0;
													for (int j = 0; j < num_foregrounds; ++ j)
														sum_pos += lyon::sqr(cbase.innerProd(integralF[j]));
													for (int j = 0; j < num_backgrounds; ++ j)
														sum_neg += lyon::sqr(cbase.innerProd(integralB[j]));
													Lvalue[ii] = sum_pos / num_foregrounds - lambda * sum_neg / num_backgrounds;
													D[ii] = 1;
													if (Lvalue[ii] > maxFunc) {
														tempL = ii;
														maxFunc = Lvalue[ii];
													}
												}

											}
										}
									}
								}
							}
					}

			}
#endif

			L.push_back(tempL);
			selected[tempL] = true;
		
		double * alpha_f = new double[num_foregrounds];
		double * alpha_b = new double[num_backgrounds];
		double * S = new double[MaxNumBases];
		CvMat ** intI = (CvMat**)malloc(sizeof(CvMat*) * MaxNumBases);
		memset(intI, 0, sizeof(CvMat*) * MaxNumBases);

		// Computing parameters
		Phi.push_back(listOfBases[L[0]].toMatrix(height, width));
		intPhi.push_back(calcIntegral(Phi[0]));
		int curtime = clock();

		double vardelta = sqrt(2. - 2. * inprod_limit);

		//#define _BOUND_PLOT
		K = 1;
		// Total Time Complexity is O(KN + K^2WH) = O(KWH(WH+K))
		while (K < MaxNumBases) {
			/**
			xmlDocPtr doc = xmlNewDoc(BAD_CAST"1.0");
			xmlNodePtr root_node = xmlNewNode(NULL, BAD_CAST"staterror");
			xmlDocSetRootElement(doc, root_node);
			static char buf[255];*/
       // Notations:
        //  * base_{k-1} = listOfBases[L[K-1]]

        // Compute alpha
        // alpha(x) = <phi_{k-1}, x> = <base_{k-1}, epsilon_{k-2}(x)>
		for (int i = 0; i < num_foregrounds; ++ i)
			alpha_f[i] = listOfBases[L[K-1]].innerProd(int_epsilon_f[i]);
        for (int i = 0; i < num_backgrounds; ++ i)
            alpha_b[i] = listOfBases[L[K-1]].innerProd(int_epsilon_b[i]);

        // u[k-1] = D[L[k-1]] = ||phi_{k-1}||^2
        double u_ksub1 = D[L[K-1]];
        double sqrt_u_ksub1 = sqrt(u_ksub1);

        // Compute S_k
        // S_k = 1/u[k-1](1/N_f\sum{alpha^2(f)}-lambda/N_b\sum{alpha^2(b)})
		double Sk = 0;
		for (int i = 0; i < num_foregrounds; ++ i)
			Sk += alpha_f[i] * alpha_f[i] / num_foregrounds / u_ksub1;
        for (int i = 0; i < num_backgrounds; ++ i)
            Sk -= lambda / num_backgrounds * alpha_b[i] * alpha_b[i] / u_ksub1;
		S[K - 1] = Sk;

        // Compute I_k
        // ita_k(x) = alpha(x)epsilon_{k-2}(x)
        // I_k = -2/sqrt(u[k-1])(1/N_f\sum{ita_k(f)}-lambda/N_b\sum{ita_k(b)})
		CvMat * Ik = cvCreateMat(height, width, NBS_CV_MAT_TYPE);
        for (int r = 0; r < Ik->rows; ++ r)
            for (int c = 0; c < Ik->cols; ++ c) {
                double curr = 0;
				for (int i = 0; i < num_foregrounds; ++ i)
					curr += alpha_f[i] * CV_MAT_ELEM(*epsilon_f[i], double, r, c) / num_foregrounds;
                for (int i = 0; i < num_backgrounds; ++ i)
                    curr -= lambda / num_backgrounds * alpha_b[i] * CV_MAT_ELEM(*epsilon_b[i], double, r, c);
                CV_MAT_ELEM(*Ik, double, r, c) = -2 * curr / sqrt_u_ksub1;
            }
        // The integral image of I_k
        CvMat * intIk = calcIntegral(Ik);
		intI[K - 1] = intIk;
		cvReleaseMat(&Ik);

            // Update epsilon_{k-2}(samples)
            // epsilon_k(x) = epsilon_{k-1}(x) - phi_k<phi_k,x>/||phi_k||^2
			for (int i = 0; i < num_foregrounds; ++ i)
				for (int r = 0; r < epsilon_f[0]->rows; ++ r)
					for (int c = 0; c < epsilon_f[0]->cols; ++ c)
						CV_MAT_ELEM(*epsilon_f[i], double, r, c)
						-= CV_MAT_ELEM(*Phi[K-1], double, r, c) * alpha_f[i] / sqrt_u_ksub1;
            for (int i = 0; i < num_backgrounds; ++ i)
                for (int r = 0; r < epsilon_b[0]->rows; ++ r)
                    for (int c = 0; c < epsilon_b[0]->cols; ++ c)
                        CV_MAT_ELEM(*epsilon_b[i], double, r, c)
                        -= CV_MAT_ELEM(*Phi[K-1], double, r, c) * alpha_b[i] / sqrt_u_ksub1;
            // Calculate the integral images of epsilon_f's and epsilon_b's.
			for (int i = 0; i < num_foregrounds; ++ i) {
				cvReleaseMat(&int_epsilon_f[i]);
				int_epsilon_f[i] = calcIntegral(epsilon_f[i]);
			}
            for (int i = 0; i < num_backgrounds; ++ i) {
                cvReleaseMat(&int_epsilon_b[i]);
                int_epsilon_b[i] = calcIntegral(epsilon_b[i]);
            }


			// Selecting the best fit Base
            maxFunc = -inf;
            tempL = -1;
            double inner_basei_phi;
            double inner_basei_Ik;
            double temp;

			for(int i = 0; i < sampN; ++ i) {
				//if (!selected[i]) {
				int index = indexOfCentBases[i];
				progress[index] = K;
			    inner_basei_phi = listOfBases[index].innerProd(intPhi[K - 1]);
                inner_basei_Ik = listOfBases[index].innerProd(intIk);
                Lvalue[index] += (inner_basei_Ik + inner_basei_phi * Sk) * inner_basei_phi;

                // Update D[i] = ||phi_i - R_Phi(phi_i)||^2
                D[index] -= inner_basei_phi * inner_basei_phi;

                temp = fabs(D[index]) < 1e-10 ? 0 : (Lvalue[index] / D[index]);

                if (temp > maxFunc) {
                    maxFunc = temp;
                    tempL = index;
                }

				sampScore[i] = temp;
			}

			if (tempL == -1) break;

			//double templ_thresh = maxFunc - fabs(maxFunc) * ratio;

			for (int ind = 0; ind < sampN; ++ ind)
				if (sampScore[ind] > maxFunc - fabs(maxFunc) * ratio) {

					HaarBase & locHaarBase = listOfBases[indexOfCentBases[ind]];
					int oW = locHaarBase.w;
					int oH = locHaarBase.h;

					double tmplthresh = oW * oH * inprod_limit * inprod_limit;

					int leftlimit, rightlimit, toplimit, botlimit;

					for (int lx = 0; lx < oW; ++ lx)
						for (int rx = lx; rx < oW; ++ rx) {
							int pw = rx - lx + 1;
							int miniH = (int)ceil(tmplthresh / pw - eps);
							for (int uy = 0; uy < oH; ++ uy)
								for (int by = uy + miniH - 1; by < oH; ++ by) {
									int ph = by - uy + 1;
									int totlimit = (int)((double)ph * ph * pw * pw / tmplthresh);
									leftlimit = (lx == 0) ? min(locHaarBase.wtl, totlimit / ph - pw) : 0;
									for (int extleft = 0; extleft <= leftlimit; ++ extleft) {
										if (rx == oW - 1) {
											rightlimit = min(imagW - locHaarBase.wtl - locHaarBase.w, totlimit / ph - pw - extleft);
										} else rightlimit = 0;
										for (int extright = 0; extright <= rightlimit; ++ extright) {
											if (uy == 0) {
												toplimit = totlimit / (pw + extleft + extright) - ph;
												toplimit = min(toplimit, locHaarBase.htl);
											} else toplimit = 0;
											for (int exttop = 0; exttop <= toplimit; ++ exttop) {
												if (by == oH - 1) {
													botlimit = (int)(totlimit / (pw + extleft + extright) - ph - exttop);
													botlimit = min(botlimit, imagH - locHaarBase.htl - locHaarBase.h);
												} else botlimit = 0;
												for (int extbot = 0; extbot <= botlimit; ++ extbot) {
													// this basis should be included in the current cluster.
													// (extleft, extright, exttop, extbot)
													int ii = calc_index_haarbase(
														locHaarBase.htl + uy - exttop,
														locHaarBase.wtl + lx - extleft,
														ph + exttop + extbot,
														pw + extleft + extright,
														imagH, imagW);

													if (progress[ii] == K || selected[ii]) continue;

													if (progress[ii] == -1) {
														HaarBase & cbase = listOfBases[ii];
														progress[ii] = 0;
														double sum_pos = 0, sum_neg = 0;
														for (int j = 0; j < num_foregrounds; ++ j)
															sum_pos += lyon::sqr(cbase.innerProd(integralF[j]));
														for (int j = 0; j < num_backgrounds; ++ j)
															sum_neg += lyon::sqr(cbase.innerProd(integralB[j]));
														Lvalue[ii] = sum_pos / num_foregrounds - lambda * sum_neg / num_backgrounds;
														D[ii] = 1;
														if (Lvalue[ii] > maxFunc) {
															tempL = ii;
															maxFunc = Lvalue[ii];
														}
													}

													// check basis-ii to update the optimal value
													for (int k = progress[ii]; k < K; ++ k) {
														HaarBase cbase = listOfBases[ii];
														inner_basei_phi = listOfBases[ii].innerProd(intPhi[k]);
														inner_basei_Ik = listOfBases[ii].innerProd(intI[k]);
														Lvalue[ii] += (inner_basei_Ik + inner_basei_phi * S[k]) * inner_basei_phi;

														// Update D[i] = ||phi_i - R_Phi(phi_i)||^2
														D[ii] -= inner_basei_phi * inner_basei_phi;

													}
													
													temp = fabs(D[ii]) < 1e-10 ? 0 : (Lvalue[ii] / D[ii]);

													if (temp > maxFunc) {
														maxFunc = temp;
														tempL = ii;
													}

													progress[ii] = K;
													
												}
											}
										}
									}
								}
						}

				}


				L.push_back(tempL);
				selected[tempL] = true;

                // Computing parameters
                CvMat * tempMat = calcPhi(L[K], Phi, intPhi);  // O(KWH)
                Phi.push_back(calcMatDiv(tempMat, sqrt(D[L[K]])));  // O(WH)
                intPhi.push_back(calcIntegral(Phi[K]));        // O(WH)
                cvReleaseMat(&tempMat);

				// Next iteration
				++ K;
		}

		lyon::settime((clock() - curtime) / (double)CLOCKS_PER_SEC);

		// Copy Answers
		this->indexOfBases.clear();
		for(int i = 0; i < (int)L.size(); ++ i)
			this->indexOfBases.push_back(listOfBases[L[i]]);
		CvMat * intTarget = calcIntegral(target);
		this->coefMatrix();
		this->coeffOfBases = compCoeff(intTarget);

		//demo();

		// Release Memories
		delete [] progress;
		delete [] D;
		delete [] selected;
		for (int i = 0; i < num_foregrounds; ++ i) {
			delete [] BF[i];
			cvReleaseMat(&int_epsilon_f[i]);
			cvReleaseMat(&epsilon_f[i]);
		}
		for (int i = 0; i < num_backgrounds; ++ i) {
			delete [] BG[i];
			cvReleaseMat(&int_epsilon_b[i]);
			cvReleaseMat(&epsilon_b[i]);
		}
		for(int i = 0; i < (int)Phi.size(); ++ i) {
			cvReleaseMat(&Phi[i]);
			cvReleaseMat(&intPhi[i]);
		}
		for (int i = 0; i < K; ++ i)
			if (intI[i])
			cvReleaseMat(&intI[i]);
		free(intI);
		delete [] S;
		delete [] alpha_f;
		delete [] alpha_b;
		delete [] Lvalue;
}

void DnbsClassGamma::computeDNBS_v1_branchbound_latest(
	int MaxNumBases,
	double lambda,
	vector<CvMat*> foreground,
	vector<CvMat*> background) {

		// Supposing frame is gray-level image with 3 channels

		// Supposing 'GenHaarBases_innerprod' completed.

		// Record the width and height
		int width = foreground[0]->width;
		int height = foreground[0]->height;

		int num_foregrounds = foreground.size();
		int num_backgrounds = background.size();

		// Definitions
		// K is the number of selected bases
		int K;
		// N is the total number of candidate bases
		int N = listOfBases.size();

		// B and D are two arrays of parameters in deciding the
		// optimized index of the next to-be-selected binary base.
		// i.e. selecting the one that maximising E_n(=|BF_n||BG_n|/D_n)
		vector<NBS_MAT_TYPE*> BF;
		for (int i = 0; i < num_foregrounds; ++ i) {
			NBS_MAT_TYPE * tempBF
				= new NBS_MAT_TYPE[N];
			BF.push_back(tempBF);
		}

		vector<NBS_MAT_TYPE*> BG;
		for (int i = 0; i < num_backgrounds; ++ i) {
			NBS_MAT_TYPE * tempBG
				= new NBS_MAT_TYPE[N];
			BG.push_back(tempBG);
		}
		NBS_MAT_TYPE * D = new NBS_MAT_TYPE[N];

		// Indicating whether the base is selected to form the subspace
		bool * selected = new bool[N];

		// Index of each selected binary bases
		vector<int> L;
		// Projection Coefficients
		vector<NBS_MAT_TYPE> Coeff;
		vector< vector<double> > Coeff_fg, Coeff_bg;

		// Phi's and their integral image
		vector<CvMat*> Phi;
		vector<CvMat*> intPhi;

		// Inner product of each Phi[K] and image F(i.e. matPatch)
		vector< vector<NBS_MAT_TYPE> > innerPhiF;
		for (int i = 0; i < num_foregrounds; ++ i) 
			innerPhiF.push_back(vector<NBS_MAT_TYPE>());
		// Inner product of each Phi[K] and image G
		vector< vector<NBS_MAT_TYPE> > innerPhiG;
		for (int i = 0; i < num_backgrounds; ++ i)
			innerPhiG.push_back(vector<NBS_MAT_TYPE>());

		// Beta's and their integral image
		vector<CvMat*> Beta;
		vector<CvMat*> intBeta;

		vector<double> fres, bres;

		// Initialization
		// Reset memories
		memset(selected, 0, sizeof(bool) * N);
		// Computing the integral image of the original image
		vector<CvMat*> integralF;
		for (int i = 0; i < num_foregrounds; ++ i)
			integralF.push_back(calcIntegral(foreground[i]));
		vector<CvMat*> integralG;
		for (int i = 0; i < num_backgrounds; ++ i)
			integralG.push_back(calcIntegral(background[i]));

		/// Representative features from the dictionary
		int sampN = indexOfCentBases.size();

		int * progress = new int[N];

		// Process of initializing parameters and get the first choice
		int tempL = 0;
		double maxFunc = -inf;
		double curFunc;
		for(int i = 0; i < N; ++ i) {
			progress[i] = 0;
			for (int j = 0; j < num_foregrounds; ++ j)
				BF[j][i] = listOfBases[i].innerProd(integralF[j]);
			for (int j = 0; j < num_backgrounds; ++ j)
				BG[j][i] = listOfBases[i].innerProd(integralG[j]);
			D[i] = 1;

			curFunc = 0;
			for (int j = 0; j < num_foregrounds; ++ j)
				curFunc += BF[j][i] * BF[j][i] / num_foregrounds;
			for (int j = 0; j < num_backgrounds; ++ j)
				curFunc -= BG[j][i] * BG[j][i] / num_backgrounds * lambda;

			if (curFunc > maxFunc) {
				tempL = i;
				maxFunc = curFunc;
			}
		}
		L.push_back(tempL);
		selected[tempL] = true;

		// Computing parameters
		Phi.push_back(listOfBases[L[0]].toMatrix(height, width));
		intPhi.push_back(calcIntegral(Phi[0]));
		Beta.push_back(cvCloneMat(Phi[0]));
		intBeta.push_back(calcIntegral(Beta[0]));
		for (int i = 0; i < num_foregrounds; ++ i) {
			Coeff_fg.push_back(vector<double>());
			Coeff_fg[i].push_back(BF[i][L[0]]);
		}
		for (int i = 0; i < num_backgrounds; ++ i) {
			Coeff_bg.push_back(vector<double>());
			Coeff_bg[i].push_back(BG[i][L[0]]);
		}

		//Coeff.push_back(BF[L[0]]);
		for (int i = 0; i < num_foregrounds; ++ i)
			innerPhiF[i].push_back(innerProd(Phi[0], foreground[i]));
		for (int i = 0; i < num_backgrounds; ++ i)
			innerPhiG[i].push_back(innerProd(Phi[0], background[i]));

		for (int i = 0; i < num_foregrounds; ++ i) {
			fres.push_back(0);
			for (int r = 0; r < foreground[i]->rows; ++ r)
				for (int c = 0; c < foreground[i]->cols; ++ c) {
					double curr = CV_MAT_ELEM(*foreground[i], double, r, c)
						- Coeff_fg[i][0] * listOfBases[L[0]].element(r, c, imagW);
					fres[i] = fres[i] + curr * curr;
				}
				fres[i] = sqrt(fres[i]);
		}

		for (int i = 0; i < num_backgrounds; ++ i) {
			bres.push_back(0);
			for (int r = 0; r < background[i]->rows; ++ r)
				for (int c = 0; c < background[i]->cols; ++ c) {
					double curr = CV_MAT_ELEM(*background[i], double, r, c)
						- Coeff_bg[i][0] * listOfBases[L[0]].element(r, c, imagW);
					bres[i] = bres[i] + curr * curr;
				}
				bres[i] = sqrt(bres[i]);
		}

		int curtime = clock();

		double * sampScore = new double[sampN];
		double * bounds = new double[sampN];

		double vardelta = sqrt(2. - 2. * inprod_limit);


		char filename[225];
		sprintf(filename, "score_stat_%.3lf.txt", inprod_limit);
		FILE * fout = fopen(filename, "w+");
		K = 1;
		// Total Time Complexity is O(KN + K^2WH) = O(KWH(WH+K))
		while (K < MaxNumBases) {
			/**
			xmlDocPtr doc = xmlNewDoc(BAD_CAST"1.0");
			xmlNodePtr root_node = xmlNewNode(NULL, BAD_CAST"staterror");
			xmlDocSetRootElement(doc, root_node);
			static char buf[255];*/

			// Selecting the best fit Base
			maxFunc = -inf;
			tempL = -1;
			double tempPhiBase;
			int optIndex;

			for(int i = 0; i < sampN; ++ i) {
				//if (!selected[i]) {
				int ii = indexOfCentBases[i];
				progress[ii] = K;
				tempPhiBase = listOfBases[ii].innerProd(intPhi[K - 1]);
				for (int j = 0; j < num_foregrounds; ++ j)
					BF[j][ii] -= innerPhiF[j][K - 1] * tempPhiBase;
				for (int j = 0; j < num_backgrounds; ++ j)
					BG[j][ii] -= innerPhiG[j][K - 1] * tempPhiBase;
				D[ii] -= tempPhiBase * tempPhiBase;

				curFunc = 0;
				for (int j = 0; j < num_foregrounds; ++ j)
					curFunc += BF[j][ii] * BF[j][ii] / num_foregrounds;
				for (int j = 0; j < num_backgrounds; ++ j)
					curFunc -= BG[j][ii] * BG[j][ii] / num_backgrounds * lambda;

				curFunc = fabs(curFunc) < sqeps ? 0 : (curFunc / D[ii]);

				sampScore[i] = curFunc;

				/**
				xmlNodePtr node = xmlNewNode(NULL, BAD_CAST "base");
				xmlAddChild(root_node, node);
				xmlNewProp(node, BAD_CAST "index", BAD_CAST itoa(i, buf, 255));
				xmlNewProp(node, BAD_CAST "func", BAD_CAST double2string(curFunc, buf, 255)); */

				if (curFunc > maxFunc) {
					maxFunc = curFunc;
					tempL = ii;
					optIndex = i;
				}

				// Calculate Bounds
				double denom = D[ii] - 2 * vardelta * sqrt(D[ii]);
				if (denom < eps) {
					bounds[i] = inf;
				} else {
					double left_term = 0;
					for (int j = 0; j < num_foregrounds; ++ j)
						left_term += lyon::sqr(BF[j][ii] + vardelta * fres[j]);
					double right_term = 0;
					for (int j = 0; j < num_backgrounds; ++ j)
						right_term += lyon::sqr(MAX(BG[j][ii] - vardelta * bres[j], 0));
					bounds[i] = left_term / num_foregrounds / denom
						- lambda * right_term / num_backgrounds / denom;
					printf("curFunc = %.6lf, bound = %.6lf\n", curFunc, bounds[i]);
				}
			}

			/**
			/// Save XML documents
			//Dumping document to stdio or file
			xmlSaveFormatFileEnc((string("exp/") + itoa(K, buf, 255) + string(".xml")).c_str(), doc, "UTF-8", 1);
			xmlFreeDoc(doc);
			xmlCleanupParser();
			xmlMemoryDump();//debug memory for regression tests */

			if (tempL == -1) break;

			int num_of_hits = 0;

			for (int ind = 0; ind < sampN; ++ ind)
				/*if (bounds[ind] > maxFunc)*/ {

					++ num_of_hits;

					double tempFunc = 0;

					for (int i = 0; i < (int)memberOfGroup[ind].size(); ++ i) {
						int ii = memberOfGroup[ind][i];

						for (int k = progress[ii]; k < K; ++ k) {
							tempPhiBase = listOfBases[ii].innerProd(intPhi[k]);
							for (int j = 0; j < num_foregrounds; ++ j)
								BF[j][ii] -= innerPhiF[j][k] * tempPhiBase;
							for (int j = 0; j < num_backgrounds; ++ j)
								BG[j][ii] -= innerPhiG[j][k] * tempPhiBase;
							D[ii] -= tempPhiBase * tempPhiBase;
						}

						progress[ii] = K;

						curFunc = 0;
						for (int j = 0; j < num_foregrounds; ++ j)
							curFunc += BF[j][ii] * BF[j][ii] / num_foregrounds;
						for (int j = 0; j < num_backgrounds; ++ j)
							curFunc -= BG[j][ii] * BG[j][ii] / num_backgrounds * lambda;

						curFunc = fabs(curFunc) < sqeps ? 0 : (curFunc / D[ii]);

						if (curFunc - sampScore[ind] > tempFunc) {
							tempFunc = curFunc - sampScore[ind];
						}

						if (curFunc > maxFunc) {
							maxFunc = curFunc;
							tempL = ii;
						}

					}

					fprintf(fout, "%d\t%.6lf\t%.6lf\t%.6lf\n", K, tempFunc, sampScore[ind], D[indexOfCentBases[ind]]);
			}


			printf("number of hits = %d\n", num_of_hits);

			L.push_back(tempL);
			selected[tempL] = true;

			// Computing parameters
			CvMat * tempMat = calcPhi(L[K], Phi, intPhi);       // O(KWH)
			Phi.push_back(calcMatDiv(tempMat, sqrt(D[L[K]])));  // O(WH)
			Beta.push_back(calcMatDiv(tempMat, D[L[K]]));       // O(WH)
			intBeta.push_back(calcIntegral(Beta[K]));           // O(WH)
			cvReleaseMat(&tempMat);
			intPhi.push_back(calcIntegral(Phi[K]));             // O(WH)
			for (int i = 0; i < num_foregrounds; ++ i)
				innerPhiF[i].push_back(innerProd(Phi[K], foreground[i]));
			for (int i = 0; i < num_backgrounds; ++ i)
				innerPhiG[i].push_back(innerProd(Phi[K], background[i]));

			// Post-processing at the end of each selection
			for (int i = 0; i < num_foregrounds; ++ i)
				Coeff_fg[i].push_back(BF[i][L[K]] / D[L[K]]);
			for (int i = 0; i < num_backgrounds; ++ i)
				Coeff_bg[i].push_back(BG[i][L[K]] / D[L[K]]);
			NBS_MAT_TYPE temp;
			for(int i = 0; i < K; ++ i) {   // O(KWH)
				temp = listOfBases[L[K]].innerProd(intBeta[i]);
				for(int r = 0; r < Beta[i]->rows; ++ r)
					for(int c = 0; c < Beta[i]->cols; ++ c)
						CV_MAT_ELEM(*Beta[i], NBS_MAT_TYPE, r, c)
						-= CV_MAT_ELEM(*Beta[K], NBS_MAT_TYPE, r, c) * temp;
				// 'temp' is real, so its conjugate is itself;
				for (int j = 0; j < num_foregrounds; ++ j)
					Coeff_fg[j][i] -= temp * Coeff_fg[j][K];
				for (int j = 0; j < num_backgrounds; ++ j)
					Coeff_bg[j][i] -= temp * Coeff_bg[j][K];
				cvReleaseMat(&intBeta[i]);
				intBeta[i] = calcIntegral(Beta[i]);
			}

			for (int i = 0; i < num_foregrounds; ++ i) {
				fres[i] = 0;
				for (int r = 0; r < foreground[i]->rows; ++ r)
					for (int c = 0; c < foreground[i]->cols; ++ c) {
						double curr = CV_MAT_ELEM(*foreground[i], double, r, c);
						for (int j = 0; j < K; ++ j) {
							curr -= Coeff_fg[i][j] * listOfBases[L[j]].element(r, c, imagW);
						}
						fres[i] = fres[i] + curr * curr;
					}
					fres[i] = sqrt(fres[i]);
			}

			for (int i = 0; i < num_backgrounds; ++ i) {
				bres[i] = 0;
				for (int r = 0; r < background[i]->rows; ++ r)
					for (int c = 0; c < background[i]->cols; ++ c) {
						double curr = CV_MAT_ELEM(*background[i], double, r, c);
						for (int j = 0; j < K; ++ j) {
							curr -= Coeff_bg[i][j] * listOfBases[L[j]].element(r, c, imagW);
						}
						bres[i] = bres[i] + curr * curr;
					}
					bres[i] = sqrt(bres[i]);
			}

			// Next iteration
			++ K;
		}

		fclose(fout);
		lyon::settime((clock() - curtime) / (double)CLOCKS_PER_SEC);

		// Copy Answers
		this->indexOfBases.clear();
		for(int i = 0; i < (int)L.size(); ++ i)
			this->indexOfBases.push_back(listOfBases[L[i]]);
		CvMat * intTarget = calcIntegral(target);
		this->coefMatrix();
		this->coeffOfBases = compCoeff(intTarget);

		//demo();

		// Release Memories
		delete [] progress;
		delete [] D;
		delete [] selected;
		for (int i = 0; i < num_foregrounds; ++ i) {
			delete [] BF[i];
		}
		for (int i = 0; i < num_backgrounds; ++ i) {
			delete [] BG[i];
		}
		for(int i = 0; i < (int)Phi.size(); ++ i) {
			cvReleaseMat(&Phi[i]);
			cvReleaseMat(&intPhi[i]);
			cvReleaseMat(&Beta[i]);
			cvReleaseMat(&intBeta[i]);
		}
}


vector<int> DnbsClassGamma::computeDNBS_v1_customFeatSet(
	int MaxNumBases,
	double lambda,
	const vector<CvMat*> & foreground,
	const vector<CvMat*> & background,
	const vector<int> & featIndex) {

		// Record the width and height
		int width = foreground[0]->width;
		int height = foreground[0]->height;

		int num_foregrounds = foreground.size();
		int num_backgrounds = background.size();

		// Definitions
		// K is the number of selected bases
		int K;
		// N is the total number of candidate bases
		int featN = featIndex.size();

		// B and D are two arrays of parameters in deciding the
		// optimized index of the next to-be-selected binary base.
		// i.e. selecting the one that maximising E_n(=|BF_n||BG_n|/D_n)
		vector<NBS_MAT_TYPE*> BF;
		for (int i = 0; i < num_foregrounds; ++ i) {
			NBS_MAT_TYPE * tempBF
				= new NBS_MAT_TYPE[featN];
			BF.push_back(tempBF);
		}

		vector<NBS_MAT_TYPE*> BG;
		for (int i = 0; i < num_backgrounds; ++ i) {
			NBS_MAT_TYPE * tempBG
				= new NBS_MAT_TYPE[featN];
			BG.push_back(tempBG);
		}
		NBS_MAT_TYPE * D = new NBS_MAT_TYPE[featN];

		// Indicating whether the base is selected to form the subspace
		bool * selected = new bool[featN];

		// Index of each selected binary bases
		vector<int> L;
		// Projection Coefficients
		vector<NBS_MAT_TYPE> Coeff;

		// Phi's and their integral image
		vector<CvMat*> Phi;
		vector<CvMat*> intPhi;

		// Inner product of each Phi[K] and image F(i.e. matPatch)
		vector< vector<NBS_MAT_TYPE> > innerPhiF;
		for (int i = 0; i < num_foregrounds; ++ i)
			innerPhiF.push_back(vector<NBS_MAT_TYPE>());
		// Inner product of each Phi[K] and image G
		vector< vector<NBS_MAT_TYPE> > innerPhiG;
		for (int i = 0; i < num_backgrounds; ++ i)
			innerPhiG.push_back(vector<NBS_MAT_TYPE>());

		// Beta's and their integral image
		vector<CvMat*> Beta;
		vector<CvMat*> intBeta;

		// Initialization
		// Reset memories
		memset(selected, 0, sizeof(bool) * featN);
		// Computing the integral image of the original image
		vector<CvMat*> integralF;
		for (int i = 0; i < num_foregrounds; ++ i)
			integralF.push_back(calcIntegral(foreground[i]));
		vector<CvMat*> integralG;
		for (int i = 0; i < num_backgrounds; ++ i)
			integralG.push_back(calcIntegral(background[i]));

		// Process of initializing parameters and get the first choice
		int tempL = 0;
		int optIndex = 0;
		double maxFunc = -inf;
		double curFunc;
		for(int i = 0; i < featN; ++ i) {
			int ii = featIndex[i];
			for (int j = 0; j < num_foregrounds; ++ j)
				BF[j][i] = listOfBases[ii].innerProd(integralF[j]);
			for (int j = 0; j < num_backgrounds; ++ j)
				BG[j][i] = listOfBases[ii].innerProd(integralG[j]);
			D[i] = 1;

			curFunc = 0;
			for (int j = 0; j < num_foregrounds; ++ j)
				curFunc += BF[j][i] * BF[j][i] / num_foregrounds;
			for (int j = 0; j < num_backgrounds; ++ j)
				curFunc -= BG[j][i] * BG[j][i] / num_backgrounds * lambda;

			if (curFunc > maxFunc) {
				tempL = ii;
				optIndex = i;
				maxFunc = curFunc;
			}
		}
		L.push_back(optIndex);
		selected[optIndex] = true;

		// Computing parameters
		Phi.push_back(listOfBases[tempL].toMatrix(height, width));
		intPhi.push_back(calcIntegral(Phi[0]));
		Beta.push_back(cvCloneMat(Phi[0]));
		intBeta.push_back(calcIntegral(Beta[0]));
		//Coeff.push_back(BF[L[0]]);
		for (int i = 0; i < num_foregrounds; ++ i)
			innerPhiF[i].push_back(innerProd(Phi[0], foreground[i]));
		for (int i = 0; i < num_backgrounds; ++ i)
			innerPhiG[i].push_back(innerProd(Phi[0], background[i]));

		int curtime = clock();

		K = 1;
		// Total Time Complexity is O(KN + K^2WH) = O(KWH(WH+K))
		while (K < MaxNumBases) {
			// Selecting the best fit Base
			maxFunc = -inf;
			tempL = -1;
			double tempPhiBase;
			for(int i = 0; i < featN; ++ i) {
				//if (!selected[i]) {
				int ii = featIndex[i];
				tempPhiBase = listOfBases[ii].innerProd(intPhi[K - 1]);
				for (int j = 0; j < num_foregrounds; ++ j)
					BF[j][i] -= innerPhiF[j][K - 1] * tempPhiBase;
				for (int j = 0; j < num_backgrounds; ++ j)
					BG[j][i] -= innerPhiG[j][K - 1] * tempPhiBase;
				D[i] -= tempPhiBase * tempPhiBase;

				curFunc = 0;
				for (int j = 0; j < num_foregrounds; ++ j)
					curFunc += BF[j][i] * BF[j][i] / num_foregrounds;
				for (int j = 0; j < num_backgrounds; ++ j)
					curFunc -= BG[j][i] * BG[j][i] / num_backgrounds * lambda;

				curFunc = fabs(curFunc) < sqeps ? 0 : (curFunc / D[i]);

				if (curFunc > maxFunc) {
					maxFunc = curFunc;
					optIndex = i;
					tempL = ii;
				}
			}

			if (tempL == -1) break;

			L.push_back(optIndex);
			selected[optIndex] = true;

			// Computing parameters
			CvMat * tempMat = calcPhi(tempL, Phi, intPhi);       // O(KWH)
			Phi.push_back(calcMatDiv(tempMat, sqrt(D[L[K]])));  // O(WH)
			Beta.push_back(calcMatDiv(tempMat, D[L[K]]));       // O(WH)
			intBeta.push_back(calcIntegral(Beta[K]));           // O(WH)
			cvReleaseMat(&tempMat);
			intPhi.push_back(calcIntegral(Phi[K]));             // O(WH)
			for (int i = 0; i < num_foregrounds; ++ i)
				innerPhiF[i].push_back(innerProd(Phi[K], foreground[i]));
			for (int i = 0; i < num_backgrounds; ++ i)
				innerPhiG[i].push_back(innerProd(Phi[K], background[i]));

			++ K;
		}

		lyon::settime((clock() - curtime) / (double)CLOCKS_PER_SEC);

		// Copy Answers
		this->indexOfBases.clear();
		for(int i = 0; i < (int)L.size(); ++ i)
			this->indexOfBases.push_back(listOfBases[featIndex[L[i]]]);
		CvMat * intTarget = calcIntegral(target);
		this->coefMatrix();
		this->coeffOfBases = compCoeff(intTarget);

		//demo();

		// Release Memories
		delete [] D;
		delete [] selected;
		for (int i = 0; i < num_foregrounds; ++ i) {
			delete [] BF[i];
		}
		for (int i = 0; i < num_backgrounds; ++ i) {
			delete [] BG[i];
		}
		for(int i = 0; i < (int)Phi.size(); ++ i) {
			cvReleaseMat(&Phi[i]);
			cvReleaseMat(&intPhi[i]);
			cvReleaseMat(&Beta[i]);
			cvReleaseMat(&intBeta[i]);
		}

		vector<int> result;
		for (int i = 0; i < (int)L.size(); ++ i)
			result.push_back(featIndex[L[i]]);
		return result;
}

vector<int> DnbsClassGamma::genDownScaleFeatures(int sampwidth, int sampheight) {
	vector<int> ret;
	for (int htl = 0; htl < imagH; htl += sampheight)
		for (int wtl = 0; wtl < imagW; wtl += sampwidth)
			for (int h = sampheight; htl + h < imagH; h += sampheight)
				for (int w = sampwidth; wtl + w < imagW; w += sampwidth) {
					ret.push_back(calc_index_haarbase(htl, wtl, h, w, imagH, imagW));
				}
				return ret;
}

#include <set>

vector<int> DnbsClassGamma::genDiffuseFeatures(const vector<int>& index, double inprod) {
	set<int> stat;
	for (int i = 0; i < (int)index.size(); ++ i) {
		HaarBase & locHaarBase = listOfBases[index[i]];
		int oW = locHaarBase.w;
		int oH = locHaarBase.h;

		double tmplthresh = oW * oH * inprod * inprod;
		//system("pause");

		int leftlimit, rightlimit, toplimit, botlimit;

		for (int lx = 0; lx < oW; ++ lx)
			for (int rx = lx; rx < oW; ++ rx) {
				int pw = rx - lx + 1;
				int miniH = (int)ceil(tmplthresh / pw - eps);
				for (int uy = 0; uy < oH; ++ uy)
					for (int by = uy + miniH - 1; by < oH; ++ by) {
						int ph = by - uy + 1;
						int totlimit = (int)((double)ph * ph * pw * pw / tmplthresh);
						leftlimit = (lx == 0) ? min(locHaarBase.wtl, totlimit / ph - pw) : 0;
						for (int extleft = 0; extleft <= leftlimit; ++ extleft) {
							if (rx == oW - 1) {
								rightlimit = min(imagW - locHaarBase.wtl - locHaarBase.w, totlimit / ph - pw - extleft);
							} else rightlimit = 0;
							for (int extright = 0; extright <= rightlimit; ++ extright) {
								if (uy == 0) {
									toplimit = totlimit / (pw + extleft + extright) - ph;
									toplimit = min(toplimit, locHaarBase.htl);
								} else toplimit = 0;
								for (int exttop = 0; exttop <= toplimit; ++ exttop) {
									if (by == oH - 1) {
										botlimit = (int)(totlimit / (pw + extleft + extright) - ph - exttop);
										botlimit = min(botlimit, imagH - locHaarBase.htl - locHaarBase.h);
									} else botlimit = 0;
									for (int extbot = 0; extbot <= botlimit; ++ extbot) {
										// this basis should be included in the current cluster.
										// (extleft, extright, exttop, extbot)
										//printf("limit: (%d, %d, %d, %d)\n", leftlimit, rightlimit, toplimit, botlimit);

										//printf("(%d, %d, %d, %d)\n", extleft, extright, exttop, extbot);
										int ii = calc_index_haarbase(
											locHaarBase.htl + uy - exttop,
											locHaarBase.wtl + lx - extleft,
											ph + exttop + extbot,
											pw + extleft + extright,
											imagH, imagW);
										stat.insert(ii);
									}
								}
							}
						}
					}
			}
	}

	vector<int> ret;
	for (set<int>::iterator it = stat.begin(); it != stat.end(); ++ it) {
		ret.push_back(*it);
	}
	printf("size of diffuse feature set = %d\n", ret.size());
	return ret;
}


void DnbsClassGamma::computeNBS_v1_cluster_v2(int MaxNumBases, double ratio) {
	CvMat * matPatch;
	// Integral image of the original image
	CvMat * intPatch;
	// K is the number of selected bases
	int K;
	// N is the total number of candidate bases
	int N = listOfBases.size();

	// B and D are two arrays of parameters in deciding the
	// optimized index of the next to-be-selected binary base.
	// i.e. selecting the one that maximising E_n(=|B_n|^2/D_n)
	NBS_MAT_TYPE * B = new NBS_MAT_TYPE[N];
	NBS_MAT_TYPE * D = new NBS_MAT_TYPE[N];

	// Indicating whether the base is selected to form the subspace
	bool * selected = new bool[N];

	// Index of each selected binary bases
	vector<int> L;
	// Projection Coefficients
	vector<NBS_MAT_TYPE> Coeff;

	// Residue[K] stands for the residue when subspace is
	// spanned by the first K binary bases.
	// i.e. Residue[K] = f - Proj_V[K](f)
	vector<NBS_MAT_TYPE> Residue;

	// Phi's and their integral image
	vector<CvMat*> Phi;
	vector<CvMat*> intPhi;

	// Inner product of each Phi[K] and image F(i.e. matPatch)
	vector<NBS_MAT_TYPE> innerPhiF;

	// Beta's and their integral image
	vector<CvMat*> Beta;
	vector<CvMat*> intBeta;

	// Temporal memories
	double curMaxE;
	int tempL;

	// Pre-processing, converting image 'patch' to type of matrix
	matPatch = cvCloneMat(target);
	// Initialization
	//printf("Initialization Begin!\n");
	double dtime = clock();
	// Reset memories
	memset(selected, false, sizeof(bool) * N);
	// Computing the integral image of the original image
	intPatch = calcIntegral(matPatch);
	// Process of initializing parameters and get the first choice

	int *progress = new int[N];
	tempL = 0;
	for(int i = 0; i < N; ++ i) {
		progress[i] = 0;
		B[i] = listOfBases[i].innerProd(intPatch);
		D[i] = 1;
		if (fabs(B[i]) > fabs(B[tempL])) tempL = i;
	}
	L.push_back(tempL);
	selected[tempL] = true;
	//printf ("Standard L[0] = %d\n", L[0]);

	// Computing parameters
	Phi.push_back(listOfBases[L[0]].toMatrix(imagH, imagW));
	intPhi.push_back(calcIntegral(Phi[0]));
	Beta.push_back(cvCloneMat(Phi[0]));
	intBeta.push_back(calcIntegral(Beta[0]));
	Coeff.push_back(B[L[0]]);
	Residue.push_back(SqrNorm(matPatch) - Sqr(Coeff[0]));
	innerPhiF.push_back(innerProd(Phi[0], matPatch));
	K = 1;

	dtime = clock();

	int sampN = indexOfCentBases.size();

	double *sampScore = new double[sampN];

	while (K < MaxNumBases && Residue[K - 1] > sigma) {
		// O(N + KWH) for each repetition
		// Selecting the best fit Base in O(N)
		curMaxE = -inf; tempL = -1;
		NBS_MAT_TYPE tempPhiBase, tempE;
		for(int i = 0; i < sampN; ++ i) {
			int ii = indexOfCentBases[i];
			progress[ii] = K;
			tempPhiBase = listOfBases[ii].innerProd(intPhi[K - 1]);
			B[ii] -= innerPhiF[K - 1] * tempPhiBase;
			D[ii] -= Sqr(tempPhiBase);
			tempE = fabs(B[i]) < eps ? 0 : (Sqr(B[i]) / D[i]);
			sampScore[i] = tempE;
			if (tempE > curMaxE) {
				curMaxE = tempE;
				tempL = i;
			}
		}

		if (tempL == -1) break;

		//double templ_thresh = maxFunc - fabs(maxFunc) * ratio;

		for (int ind = 0; ind < sampN; ++ ind)
			if (sampScore[ind] > curMaxE - fabs(curMaxE) * ratio) {

				HaarBase & locHaarBase = listOfBases[indexOfCentBases[ind]];
				int oW = locHaarBase.w;
				int oH = locHaarBase.h;

				double tmplthresh = oW * oH * inprod_limit * inprod_limit;

				int leftlimit, rightlimit, toplimit, botlimit;

				for (int lx = 0; lx < oW; ++ lx)
					for (int rx = lx; rx < oW; ++ rx) {
						int pw = rx - lx + 1;
						int miniH = (int)ceil(tmplthresh / pw - eps);
						for (int uy = 0; uy < oH; ++ uy)
							for (int by = uy + miniH - 1; by < oH; ++ by) {
								int ph = by - uy + 1;
								int totlimit = (int)((double)ph * ph * pw * pw / tmplthresh);
								leftlimit = (lx == 0) ? min(locHaarBase.wtl, totlimit / ph - pw) : 0;
								for (int extleft = 0; extleft <= leftlimit; ++ extleft) {
									if (rx == oW - 1) {
										rightlimit = min(imagW - locHaarBase.wtl - locHaarBase.w, totlimit / ph - pw - extleft);
									} else rightlimit = 0;
									for (int extright = 0; extright <= rightlimit; ++ extright) {
										if (uy == 0) {
											toplimit = totlimit / (pw + extleft + extright) - ph;
											toplimit = min(toplimit, locHaarBase.htl);
										} else toplimit = 0;
										for (int exttop = 0; exttop <= toplimit; ++ exttop) {
											if (by == oH - 1) {
												botlimit = (int)(totlimit / (pw + extleft + extright) - ph - exttop);
												botlimit = min(botlimit, imagH - locHaarBase.htl - locHaarBase.h);
											} else botlimit = 0;
											for (int extbot = 0; extbot <= botlimit; ++ extbot) {
												// this basis should be included in the current cluster.
												// (extleft, extright, exttop, extbot)
												int ii = calc_index_haarbase(
													locHaarBase.htl + uy - exttop,
													locHaarBase.wtl + lx - extleft,
													ph + exttop + extbot,
													pw + extleft + extright,
													imagH, imagW);

												if (progress[ii] == K) continue;

												// check basis-ii to update the optimal value
												for (int k = progress[ii]; k < K; ++ k) {
													tempPhiBase = listOfBases[ii].innerProd(intPhi[k]);
													B[ii] -= innerPhiF[k] * tempPhiBase;
													D[ii] -= tempPhiBase * tempPhiBase;
												}

												progress[ii] = K;

												tempE = fabs(B[ii]) < eps ? 0 : (Sqr(B[ii]) / D[ii]);
												if (tempE > curMaxE) {
													curMaxE = tempE;
													tempL = ii;
												}
											}
										}
									}
								}
							}
					}

			}

		L.push_back(tempL);
		selected[tempL] = true;

		// Computing parameters
		Residue.push_back(Residue[K - 1] - curMaxE);        // O(1)
		CvMat * tempMat = calcPhi(L[K], Phi, intPhi);       // O(KWH)
		Phi.push_back(calcMatDiv(tempMat, sqrt(D[L[K]])));  // O(WH)
		Beta.push_back(calcMatDiv(tempMat, D[L[K]]));       // O(WH)
		intBeta.push_back(calcIntegral(Beta[K]));           // O(WH)
		cvReleaseMat(&tempMat);
		intPhi.push_back(calcIntegral(Phi[K]));             // O(WH)
		innerPhiF.push_back(innerProd(Phi[K], matPatch));   // O(WH)
		Coeff.push_back(B[L[K]] / D[L[K]]);                 // O(1)
		//printf("Rep.%d Residue[%d] = %.10lf, Coeff[%d] = %.3lf\n", K, K, sqrt(Residue[K]), K, Coeff[K]);

		// Post-processing at the end of each selection
		NBS_MAT_TYPE temp;
		for(int i = 0; i < K; ++ i) {   // O(KWH)
			temp = listOfBases[L[K]].innerProd(intBeta[i]);
			for(int r = 0; r < Beta[i]->rows; ++ r)
				for(int c = 0; c < Beta[i]->cols; ++ c)
					CV_MAT_ELEM(*Beta[i], NBS_MAT_TYPE, r, c)
					-= CV_MAT_ELEM(*Beta[K], NBS_MAT_TYPE, r, c) * temp;
			// 'temp' is real, so its conjugate is itself;
			Coeff[i] -= temp * Coeff[K];
			cvReleaseMat(&intBeta[i]);
			intBeta[i] = calcIntegral(Beta[i]);
		}

		//fprintf(outfile, "%d\t%.3lf\t%.3lf\n", K + 1, sqrt(Residue[K] / imagW / imagH), (clock() - dtime) / (double)CLOCKS_PER_SEC);
		// Now, get the next binary base
		++ K;
	}

	// Copy Answers
	this->indexOfBases.clear();
	for(int i = 0; i < (int)L.size(); ++ i)
		this->indexOfBases.push_back(listOfBases[L[i]]);
	this->coeffOfBases = Coeff;

	// Release Memories
	delete [] progress;
	delete [] sampScore;
	cvReleaseMat(&matPatch);
	cvReleaseMat(&intPatch);
	delete [] B;
	delete [] D;
	delete [] selected;
	for(int i = 0; i < (int)Phi.size(); ++ i) {
		cvReleaseMat(&Phi[i]);
		cvReleaseMat(&intPhi[i]);
		cvReleaseMat(&Beta[i]);
		cvReleaseMat(&intBeta[i]);
	}
}

void DnbsClassGamma::sample_backgrounds_supress(CvRect rect,
											   int MaxNumBackgrounds,
											   int distance_threshold,
											   int x_min, int x_max,
											   int y_min, int y_max,
											   vector< CvRect > & pnt,
											   vector< pair<double, int> > & arr,
											   double ** weight,
											   vector < CvRect > & samples,
											   double nonoverlap_ratio) 
{
	// Iteratively sample 'MaxNumBackgrounds' backgrounds
	// (1) by avoiding the similarity (SSD/NCC) of any pair
	// less than a given threshold.
	// or (2) by avoiding the distance between any pair
	// of samples' locations less than a give threshold.
	// (Here we choose (2))

	samples.clear();

	// Clear around the foreground
	for (int w = (int)ceil(nonoverlap_ratio * rect.width); w <= rect.width; ++ w)
		for (int h = (int)ceil(nonoverlap_ratio * rect.width * rect.height / w); h <= rect.height; ++ h) {
			int x, y;
			x = rect.x + w - rect.width;
			if (x >= x_min && x < x_max) {
				y = rect.y + h - rect.height;
				if (y >= y_min && y < y_max)
					weight[x - x_min][y - y_min] = -1;
				y = rect.y + rect.height - h;
				if (y < y_max && y >= y_min)
					weight[x - x_min][y - y_min] = -1;
			}
			x = rect.x + rect.width - w;
			if (x < x_max && x >= x_min) {
				y = rect.y + h - rect.height;
				if (y >= y_min && y < y_max)
					weight[x - x_min][y - y_min] = -1;
				y = rect.y + rect.height - h;
				if (y < y_max && y >= y_min)
					weight[x - x_min][y - y_min] = -1;
			}
		}

	for (int i = 0; i < (int)arr.size() && (int)samples.size() < MaxNumBackgrounds; ++ i) {
		if (weight[pnt[arr[i].second].x - x_min][pnt[arr[i].second].y - y_min] != -1)
			samples.push_back(pnt[arr[i].second]);
		//else continue;
		int dest_x = pnt[arr[i].second].x;
		int dest_y = pnt[arr[i].second].y;
		for (int w = (int)ceil(nonoverlap_ratio * rect.width); w <= rect.width; ++ w)
			for (int h = (int)ceil(nonoverlap_ratio * rect.width * rect.height / w); h <= rect.height; ++ h) {
				int x, y;
				x = dest_x + w - rect.width;
				if (x >= x_min && x < x_max) {
					y = dest_y + h - rect.height;
					if (y >= y_min && y < y_max)
						weight[x - x_min][y - y_min] = -1;
					y = dest_y + rect.height - h;
					if (y < y_max && y >= y_min)
						weight[x - x_min][y - y_min] = -1;
				}
				x = dest_x + rect.width - w;
				if (x < x_max && x >= x_min) {
					y = dest_y + h - rect.height;
					if (y >= y_min && y < y_max)
						weight[x - x_min][y - y_min] = -1;
					y = dest_y + rect.height - h;
					if (y < y_max && y >= y_min)
						weight[x - x_min][y - y_min] = -1;
				}
			}
	}
}
