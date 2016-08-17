/*
 * funneltest.h
 *
 *  Created on: Aug 17, 2016
 *      Author: elias
 */

#ifndef TESTS_FUNNELTEST_H_
#define TESTS_FUNNELTEST_H_

int funnelBasicTest() {
	Coords ca(0.0, 0.0);
	Coords cl(0.0, 1.0);
	Coords cr(1.0, 0.0);
	Funnel fun(&cl, &ca, &cr);

	Edge base = fun.getBase();
	Coords* apex = fun.getApex();

	if (apex != &ca) {
		std::cout << "apex: " << apex << std::endl;
		std::cout << "orig: " << &ca << std::endl;
		std::cout << "apex was not at same address as original apex"
				<< std::endl;
		return 1;
	}

	if (base.getL() != &cl) {
		return 2;
	}
	if (base.getR() != &cr) {
		return 3;
	}

	return 0;

}

int funnelBaseTest() {
	Coords ca(0.0, 0.0);
	Coords cl1(0.5, 1.0);
	Coords cl2(0.4, 2.0);
	Coords cr1(1.0, 0.5);
	Coords cr2(2.0, 0.4);

	Funnel fun(&cl1, &ca, &cr1);
	fun.addL(&cl2);
	if (fun.getBase().getL() != &cl2 or fun.getBase().getR() != &cr1) {
		return 1;
	}
	fun.addR(&cr2);
	if (fun.getBase().getL() != &cl2 or fun.getBase().getR() != &cr2) {
		return 2;
	}
	return 0;

}

int testFunnel() {
	if (funnelBasicTest() != 0) {
		return funnelBasicTest();
	}
	if (funnelBaseTest() != 0) {
		return funnelBaseTest();
	}
	return 0;
}

#endif /* TESTS_FUNNELTEST_H_ */
