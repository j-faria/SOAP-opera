#include "Data.h"
#include <iostream>
#include <fstream>
#include <string>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

#define MAXBUFSIZE  ((int) 1e6)
typedef vector <double> record_t;
typedef vector <record_t> data_t;

Data Data::instance;

Data::Data(){}

MatrixXd readMatrix(const char *filename)
{
    int cols = 0, rows = 0;
    double buff[MAXBUFSIZE];

    // Read numbers from file into buffer.
    ifstream infile;
    infile.open(filename);
    while (!infile.eof())
        {
        string line;
        getline(infile, line);

        int temp_cols = 0;
        stringstream stream(line);
        while(! stream.eof())
            stream >> buff[cols*rows+temp_cols++];

        if (temp_cols == 0)
            continue;

        if (cols == 0)
            cols = temp_cols;

        rows++;
        }

    infile.close();

    rows--;

    cout << "gello!" << endl;

    // Populate matrix with numbers.
    MatrixXd result(rows,cols);
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            result(i,j) = buff[ cols*i+j ];

    return result;
};



//-----------------------------------------------------------------------------
// Let's overload the stream input operator to read a list of CSV fields (which a CSV record).
istream& operator >> ( istream& ins, record_t& record )
  {
  // make sure that the returned record contains only the stuff we read now
  record.clear();

  // read the entire line into a string (a CSV record is terminated by a newline)
  string line;
  getline( ins, line );

  // now we'll use a stringstream to separate the fields out of the line
  stringstream ss( line );

  // convert each field to a double and 
  // add the newly-converted field to the end of the record
  double f;
  while (ss >> f)
    record.push_back(f);
  
  // Now we have read a single line, converted into a list of fields, converted the fields
  // from strings to doubles, and stored the results in the argument record, so
  // we just return the argument stream as required for this kind of input overload function.
  return ins;
  }

//-----------------------------------------------------------------------------
// Let's likewise overload the stream input operator to read a list of CSV records.
// This time it is a little easier, just because we only need to worry about reading
// records, and not fields.
istream& operator >> ( istream& ins, data_t& data )
  {
  // used to hold the header
  data_t header;
  // make sure that the returned data only contains the CSV data we read here
  data.clear();

  // For every record we can read from the file, append it to our resulting data
  // except if it's in the header
  int i = 0;
  record_t record;
  while (ins >> record)
    {
      //if (i==0 || i==1) header.push_back(record);
      data.push_back(record);
      //  i++;
    }

  // Again, return the argument stream
  return ins;  
  }


void Data::load(const char* filename, const char* units, int skip)
  /* 
  Read in tab/space separated file `filename` with columns
  time  vrad  error
  ...   ...   ...
  where vrad and error are in `units` (either "kms" or "ms")
  */
  {

  data_t data;
  ccfdata = readMatrix("test_data.txt");


  // Empty the vectors
  t.clear();
  y.clear();
  sig.clear();

  // Read the file into the data container
  ifstream infile( filename );
  infile >> data;
  //operator>>(infile, data, skip);

  // Complain if something went wrong.
  if (!infile.eof())
  {
    printf("Could not read data file (%s)!\n", filename);
    exit(1);
  }

  infile.close();

  double factor = 1.;
  if(units == "kms") factor = 1E3;

  for (unsigned n = 0; n < data.size(); n++)
    {
      if (n<skip) continue;
      t.push_back(data[n][0]);
      y.push_back(data[n][1] * factor);
      sig.push_back(data[n][2] * factor);
    }

  // How many points did we read?
  printf("# Loaded %d data points from file %s\n", t.size(), filename);

  for(unsigned i=0; i<data.size(); i++)
  {
      if (t[i] > 57170.)
      {
          index_fibers = i;
          break;
      }
  }

  }