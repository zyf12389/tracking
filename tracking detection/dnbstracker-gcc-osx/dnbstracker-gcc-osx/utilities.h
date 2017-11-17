#ifndef _H_UTILITIES_
#define _H_UTILITIES_
#include <stdexcept>
#include <vector>
#include <string>
#include <sys/stat.h>
#include <dirent.h>

const float mm_per_pix = 0.1042;
const float ref_dist = 120;
const float factor = ref_dist / 2 / mm_per_pix;

void convert_from_image_to_world(float& x, float& y, float z);
void convert_from_world_to_image(float& x, float& y, float z);

#define sqr(x) ((x)*(x))
float randomf();
int random(int N);
int str2int(const std::string& s);
std::string num2str(int n);
std::vector<int> sample_K_from_N(int K, int N);
std::pair<int, int> convert_str2ints(const std::string& str, const std::string& separate);

std::string get_suffix(std::string s);
std::string change_ext(std::string s, std::string ext);
std::string remove_ext(std::string s);
bool create_dir(std::string folder);
void remove_dir(std::string folder);
void create_parent_dir(std::string filename);
void touch_file(std::string filename);
void touch_file(std::string filename, std::string message);
void remove_file(std::string filename);
bool file_exists(std::string filename);
bool dir_exists(std::string filename);

// file reading

void read_linewise_file(std::string filename, std::vector<std::string>& output);

// data sampling
void sample_depth_dataset(std::string dir_depth, std::string dir_sample, const std::vector<std::string>& labels, int num_samples);
void sample_depth_file(std::string file_depth, std::string file_sample, int num_samples);
#endif
