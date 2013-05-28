#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

std::string pts_filename(const char* filename)
{
  std::string path(filename);
  std::size_t found = path.find_last_of("/\\");
  std::string file = path.substr(found+1, path.length()-4);
  file.append(".pts");
  return file;
}

void off_to_pts(const char* filename/*.off file*/)
{
  std::string line;
  std::ifstream ifs(filename);
  std::ofstream ofs(pts_filename(filename).data());
  if(ifs.is_open())
  {
    bool nbv_found = false;
    std::size_t nbv = 0;
    int i = 0;
    while(ifs.good())
    {
      std::getline(ifs,line);
      if(line.find("OFF") != std::string::npos)
        continue;
      else if(line.find("#") != std::string::npos)//comments
        continue;
      else if(!nbv_found)
      {
        sscanf(line.c_str(),"%d",&nbv);
        std::cout << "Nb poles is : " << nbv << std::endl;
        nbv_found = true;
      }
      else if(i++ < nbv)
      {
        double x,y,z;
        std::istringstream liness(line);
        liness >> x >> y >> z;
        ofs << x << " " << y << " "<< z << std::endl;
      }
      else
        break;
    }
    ifs.close();
    ofs.close();
  }
  else 
    std::cout << "Unable to open the file : " << filename << "\n"; 
}



int main(int argc, char *argv[])
{
  if(argc != 2)
  {
    std::cout << "No file found" << std::endl;
    return -1;
  }
  
  off_to_pts(argv[1]);

  return 0;
}