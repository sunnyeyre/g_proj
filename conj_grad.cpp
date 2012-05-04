//
//  conj_grad.cpp
//  
//
//  Created by Sunling Yang on 1/25/12.
//

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <cassert>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cmath>

#include <iostream>
#include "linalg.h"
#include "gsl/gsl_multimin.h"
#include "gsl/gsl_errno.h"
#include "gsl/gsl_math.h"
#include "gsl/gsl_min.h"

#define PI 3.14159265
#define g 9.80665

using namespace std;

