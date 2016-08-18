/*
 * test.h
 *
 *  Created on: Aug 17, 2016
 *      Author: elias
 */

#ifndef TEST_H_
#define TEST_H_

#include <iostream>
#include "../Edge.h"
#include "../Funnel.h"
#include "funneltest.h"
#include "basictest.h"



void test() {
	int t = 0;
	std::cout << "!!!Starting tests!!!" << std::endl; // prints !!!Hello World!!!
	if (testCoords() != 0) {
		std::cout << "testCoords failed" << std::endl;
	}
	if (testEdge() != 0) {
		std::cout << "testEdge failed" << std::endl;
	}
	if (testDist() != 0) {
		std::cout << "testDist failed" << std::endl;
	}
	t = testIsRight();
	if (t != 0) {
		std::cout << "testIsRight failed with value: " << t << std::endl;
	}
	if (testFunnel() != 0) {
		std::cout << "testFunnel failed" << std::endl;
	}

	std::cout << "!!!Tests finished!!!" << std::endl; // prints !!!Hello World!!!

}

#endif /* TEST_H_ */
