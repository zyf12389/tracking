#include "utilities.h"
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <unistd.h>

int str2int(const std::string& s) {
	int val;
	sscanf(s.c_str(), "%d", &val);
	return val;
}

void convert_from_image_to_world(float& x, float& y, float z) {
	x = x * z / factor;
	y = y * z / factor;
}

void convert_from_world_to_image(float& x, float& y, float z) {
	x = x / z * factor;
	y = y / z * factor;
}

std::pair<int, int> convert_str2ints(const std::string& str, const std::string& separate) {
	int index = (int)str.find(separate);
	return std::make_pair(str2int(str.substr(0, index)), str2int(str.substr(index+1)));
}

bool create_dir(std::string folder) {
	struct stat st;
	if (stat(folder.c_str(), &st) != 0 || (st.st_mode & S_IFDIR) == 0)
		mkdir(folder.c_str(), S_IRWXG | S_IRWXU);
	else return false;
	return true;
}

void remove_dir(std::string folder) {
	rmdir(folder.c_str());
}

void remove_file(std::string file) {
	if (!file_exists(file))
		return;
	if (remove(file.c_str()) != 0)
		throw std::runtime_error("Unable to remove file");
}

std::string num2str(int n) {
	std::stringstream ss;
	ss << n;
	return ss.str();
}

float randomf() {
	return (float)(rand() / ((double)RAND_MAX));
}

int random(int N) {
	return (int)((rand() / ((double)RAND_MAX + 1)) * N);
}

std::vector<int> sample_K_from_N(int K, int N) {
	std::vector<int> sample;
	for (int i = 0; i < K; ++ i) {
		int offset = random(N - i);
		int j = 0;
		while (j < i && offset >= sample[j]) {
			++ offset;
			++ j;
		}
		while (j < i) {
			std::swap(offset, sample[j]);
			++ j;
		}
		sample.push_back(offset);
	}
	return sample;
}

std::string get_suffix(std::string s) {
	int index = (int)s.find_last_of(".");
	return s.substr(index);
}

std::string change_ext(std::string s, std::string ext) {
	int index = (int)s.find_last_of(".");
	if (ext[0] != '.') ext = "." + ext;
	return s.substr(0, index) + ext;
}

std::string remove_ext(std::string s) {
	int index = (int)s.find_last_of(".");
	return s.substr(0, index);
}
void create_parent_dir(std::string filename) {
	int index = (int)filename.find_last_of("/");
	std::string folder = filename.substr(0, index);
	fprintf(stderr, "folder: %s\n", folder.c_str());
	create_dir(folder);
}

void touch_file(std::string filename) {
	FILE* f = fopen(filename.c_str(), "w+");
	fclose(f);
}

void touch_file(std::string filename, std::string message) {
	FILE* f = fopen(filename.c_str(), "w+");
	fprintf(f, "%s\n", message.c_str());
	fclose(f);
}
bool file_exists(std::string filename) {
	struct stat st;
	return (stat(filename.c_str(), &st) == 0) && (st.st_mode & S_IFREG);
}

bool dir_exists(std::string folder) {
	struct stat st;
	if (stat(folder.c_str(), &st) != 0 || (st.st_mode & S_IFDIR) == 0)
		return false;
	return true;
}

void sample_depth_file(std::string file_depth, std::string file_sample, int num_samples) {
	// Read depth image
	FILE* fin = fopen(file_depth.c_str(), "rb");
	int nrows, ncols;
	fread(&nrows, sizeof(int), 1, fin);
	fread(&ncols, sizeof(int), 1, fin);
	fclose(fin);
	fprintf(stderr, "nrows = %d, ncols = %d\n", nrows, ncols);
	std::vector<int> samples = sample_K_from_N(num_samples, nrows * ncols);
	// each line of sample file contains row number and column number for each sample
	FILE* fout = fopen(file_sample.c_str(), "w+");
	if (fout == NULL) {
		throw std::runtime_error("Unable to open sample file " + file_sample);
	}

	for (std::vector<int>::iterator it = samples.begin(); it != samples.end(); ++ it)
		fprintf(fout, "%d %d\n", *it / ncols, *it % ncols);
	fclose(fout);
}

void sample_depth_dataset(std::string dir_depth, std::string dir_sample, const std::vector<std::string>& labels, int num_samples) {
	if (file_exists(dir_sample + "/.complete")) {
		fprintf(stderr, "Data have been sampled already, no need to resample.\n");
		return;
	}
	create_dir(dir_sample);
	std::string file_s, file_out;
	for (std::vector<std::string>::const_iterator it = labels.begin(); it != labels.end(); ++ it) {
		file_s = dir_depth + "/" + *it + ".depth";
		create_parent_dir(dir_sample + "/" + *it);
		file_out = dir_sample + "/" + *it + ".txt";
		sample_depth_file(file_s, file_out, num_samples);
	}
	touch_file(dir_sample + "/.complete", num2str(num_samples));
	fprintf(stderr, "Data sample finished.\n");
}

void read_linewise_file(std::string filename, std::vector<std::string>& output) {
	std::string elem;
	std::ifstream fs;
	output.clear();
	try {
		fs.open(filename.c_str());
		while (fs >> elem)
			output.push_back(elem);
		fs.close();
	} catch(std::ifstream::failure e) {
		throw std::runtime_error("Unable to open file");
	}
	std::vector<std::string>(output).swap(output);
}
