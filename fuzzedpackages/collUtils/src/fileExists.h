/*
 * fileExists.h
 *
 *  Created on: Jul 26, 2014
 *      Author: kaiyin
 */

#pragma once

#include <cstdlib>
#include <unistd.h>
#include <string>
#include <cstdio>
#include <iostream>
#include <sstream>




void fileExists (const std::string& name);
void fileExists(const std::string& name, bool exitIfNotExist, bool exitIfExist);
bool fileExists(const std::string& name, bool rmIfExist);
