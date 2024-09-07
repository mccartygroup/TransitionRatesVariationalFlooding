/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2023 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "bias/Bias.h"
#include "bias/ActionRegister.h"
#include "core/ActionSet.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include "tools/Grid.h"
#include "tools/OFile.h"
#include "tools/IFile.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <vector>
#include <limits>
#include <string>
#include <algorithm>
#include "spline.cpp"

using namespace std;

namespace PLMD {
namespace bias {

//+PLUMEDOC BIAS LOGFLOODING
/*

*/
//+ENDPLUMEDOC

class LogFlooding : public Bias {
private:
    unsigned int ncv_;
    double temp_;
    double kB_temp_;
    double beta_;

    Spline s;
    int _time_step_;
    int iteration;
    int n_data;

    bool _write_flooding_;
    double _a_;
    double _lambda_;

    string _vesbias_filename_;
    string _wbias_filename_;
    int  _wbias_stride_;
    double _wbias_grid_max_;
    double _wbias_grid_min_;
    double _wbias_grid_bin_;
    double _wbias_grid_size_;

    vector<double> cv_data;
    vector<double> ves_data;

    double x, y;

    string line;
    vector<double> aa;
    vector<double> bb;
    vector<double> cc;


    double get_max_value(vector<double>);
    double get_FF_value(double, double);
    double get_fermi_value(double, double, double);
    double get_fermi_derivative(double, double, double, double);


public:
  explicit LogFlooding(const ActionOptions&);
  ~LogFlooding();
  void calculate();
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(LogFlooding,"LOGFLOODING")

void LogFlooding::registerKeywords(Keywords& keys) {
  Bias::registerKeywords(keys);
  ActionWithValue::useCustomisableComponents(keys);
  keys.use("ARG");
  // Should be _bias below
//  keys.addOutputComponent("bias","default","the instantaneous value of the bias");
  keys.addOutputComponent("swbias","default","value of the switching function");
  keys.addOutputComponent("force","default","value of the force");

  keys.add("compulsory","TEMP","Temperature.");
  keys.add("compulsory","VES_BIAS_FILE","ves bias file to be read");
  keys.add("compulsory" ,"A_VALUE","value of a in the time");
  keys.add("optional","FERMI_LAMBDA","The value of lambda in the switching function. default is 100");
  keys.add("optional","OUTPUT_BIAS_FILE","File name to output the bias");
  keys.add("optional","OUTPUT_BIAS_FILE_STRIDE","Stride to output bias file");
  keys.add("optional","OUTPUT_BIAS_GRID_MIN","grid min for output bias. default is -1.0");
  keys.add("optional","OUTPUT_BIAS_GRID_MAX","grid max for output bias. default is 1.0");
  keys.add("optional","OUTPUT_BIAS_GRID_BIN","output bias grid bin. default is 100.0");
}

LogFlooding::~LogFlooding(){

}

LogFlooding::LogFlooding(const ActionOptions&ao):
PLUMED_BIAS_INIT(ao)
{
  ncv_=getNumberOfArguments();
  // To do: can only use one CV
  // Read temperature
  parse("TEMP",temp_);
  beta_=1.0/(plumed.getAtoms().getKBoltzmann()*temp_);
  kB_temp_ = plumed.getAtoms().getKBoltzmann()*temp_;
  log.printf("Temperature set %f %f %f\n",temp_, beta_, kB_temp_);
  // Read VES bias file name
  _vesbias_filename_="";
  parse("VES_BIAS_FILE",_vesbias_filename_);
  // For pulling the bias
  parse("A_VALUE",_a_);
  log.printf("Value of a in the switching function is set to %f\n",_a_);
  _lambda_=100;
  parse("FERMI_LAMBDA",_lambda_);
  log.printf("Value of lambda in the switching function is set to %f\n",_lambda_);

  _write_flooding_ = false;
  _wbias_filename_="";
  parse("OUTPUT_BIAS_FILE",_wbias_filename_);
  if(_wbias_filename_.length() > 0){_write_flooding_=true;}

  _wbias_stride_= 1;
  parse("OUTPUT_BIAS_FILE_STRIDE",_wbias_stride_);
  _wbias_grid_min_=-1.0;
  parse("OUTPUT_BIAS_GRID_MIN",_wbias_grid_min_);
  _wbias_grid_max_=1.0;
  parse("OUTPUT_BIAS_GRID_MAX",_wbias_grid_max_);
  _wbias_grid_bin_=100.0;
  parse("OUTPUT_BIAS_GRID_BIN",_wbias_grid_bin_);
  // Everything has to be parsed by now.
  checkRead();
  // add one bias for each argument
  //addComponent("bias"); componentIsNotPeriodic("bias");
  //valueBias=getPntrToComponent("bias");
  addComponent("swbias"); componentIsNotPeriodic("swbias");
  //valueSwBias=getPntrToComponent("swbias");
  addComponent("force"); componentIsNotPeriodic("force");
  //valueForceTot=getPntrToComponent("force");

  // counter
  _time_step_ = 1;
  iteration = 0;

  // grid bin size:
  _wbias_grid_size_ = (_wbias_grid_max_ - _wbias_grid_min_)/_wbias_grid_bin_;

  // Now read VES bias:
  log.printf("Reading VES bias file\n");
  ifstream data_file(_vesbias_filename_);
  if (!data_file) {error("Unable to read VES bias file\n");}

  while(getline(data_file, line)){
    if(line.find("#!") == 0){
      continue;
    }
  istringstream iss(line);
  if (iss >> x >> y) {
    cv_data.push_back(x);
    ves_data.push_back(y);
   }
  }
  data_file.close();

  // Fit VES bias data file to spline:

  n_data = cv_data.size();
  aa.resize(n_data);
  bb.resize(n_data);
  cc.resize(n_data);

  aa.assign(n_data,0.0);
  bb.assign(n_data,0.0);
  cc.assign(n_data,0.0);

  s.set_spline(cv_data,ves_data, aa, bb,cc);
  log.printf("value of n_data %f\n",n_data);
  //log.printf("aa bb cc %f %f %f\n",aa, bb, cc);
  // TO DO: add restart option here
}

// getting value of maximum bias
double LogFlooding::get_max_value(vector<double> y_data) {
  auto max_iter = std::max_element(y_data.begin(), y_data.end());

  double Vmax = *max_iter;
  return Vmax;
}

// estimate F(s) from -V(s)
double LogFlooding::get_FF_value(double vs, double Vmax){
  return -1.0*vs + Vmax;
}

// get switching function
double LogFlooding::get_fermi_value(double Fs, double a, double lambda){
  return 1.0 / (1.0 + exp(lambda*(Fs - a)));
}

// get derivative of switching function
double LogFlooding::get_fermi_derivative(double Fs, double a, double lambda, double dydx){
  double f1;
  double f2;
  f1 = exp(lambda*(Fs-a));
  f2 = 1.0/((1.0+f1)*(1.0+f1));
  return lambda*f1*f2*dydx;
}


void LogFlooding::calculate() {
  double Vmax;
  double FF_prime;
  double SW_value;
  double bias;
  double force_;
  double ves_value;
  double dx;

  // Get value of the CV (Only one CV allowed at the moment)
  double cv_value;
  cv_value=getArgument(0);

  Vmax = get_max_value(ves_data);

  // get time
  double t;
  t = 120.0*log10(1 + _time_step_ * _a_);

  // the value of the VES bias
  ves_value = s.seval(cv_data,ves_data,aa,bb,cc,cv_value);

  // The value of the free energy
  FF_prime = get_FF_value(ves_value,Vmax);

  // The new switching function value
  SW_value = get_fermi_value(FF_prime,t,_lambda_);
  getPntrToComponent("swbias")->set(SW_value);
  // The value of the bias after applying the new switching function
  bias = (-FF_prime + t)*SW_value;
  getPntrToComponent("bias")->set(bias);

  // The force
  dx = s.seval_derivative(cv_data,ves_data,aa,bb,cc,cv_value);

  force_ = -SW_value*dx - (-FF_prime + t)*get_fermi_derivative(FF_prime,t,_lambda_,dx);
  getPntrToComponent("force")->set(force_);
  for(unsigned i=0;i< getNumberOfArguments() ;++i){
//  log<<"BIAS "<<val<<"\n";
  setOutputForce(i,force_);
  }

  _time_step_ += 1;
  if(_write_flooding_){
    if(_time_step_%_wbias_stride_==0){
    string bias_iter = to_string(iteration);
    ofstream outFile;
    string bias_name = _wbias_filename_ + bias_iter;
    cout << iteration << " " << bias_iter  <<  " " << bias_name <<  endl;
    outFile.open(bias_name); // open file for writing

    for (double x = _wbias_grid_min_; x <= _wbias_grid_max_; x += _wbias_grid_size_){
      ves_value = s.seval(cv_data,ves_data,aa,bb,cc,x);
      FF_prime = get_FF_value(ves_value,Vmax);
      SW_value = get_fermi_value(FF_prime,t,_lambda_);
      double y = (-FF_prime + t)*SW_value;
      // The force
      dx = s.seval_derivative(cv_data,ves_data,aa,bb,cc,x);

      double dy = -SW_value*dx - (-FF_prime + t)*get_fermi_derivative(FF_prime,t,_lambda_,dx);
      outFile << x << " " << y << " " << dy <<  endl;
    }
    outFile.close();
    iteration = iteration + 1;
  }
  }
  }
}
}
