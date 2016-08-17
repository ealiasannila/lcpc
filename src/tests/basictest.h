/*
 * basictest.h
 *
 *  Created on: Aug 17, 2016
 *      Author: elias
 */

#ifndef TESTS_BASICTEST_H_
#define TESTS_BASICTEST_H_

int testCoords() {
	Coords c(2.0, 1.0);
	if (c.getX() == 2.0 and c.getY() == 1.0) {
		return 0;
	}
	return 1;
}

int testEdge() {
	Coords c1(0.0, 0.0);
	Coords c2(2.0, 2.0);
	Edge e(&c1, &c2);
	Coords* cr = e.getR();
	Coords* cl = e.getL();
	if (cl != &c1) {
		return 1;
	}
	if (cr != &c2) {
		return 2;
	}
	return 0;
}

int testDist() {
	Coords c1(0.0, 0.0);
	Coords c2(2.0, 2.0);
	Coords c3(0.0, 1.0);

	if (eucDist(c1, c3) != 1 or eucDist(c3, c1) != 1) {
		return 1;
	}
	if (eucDist(c1, c2) != sqrt(8) or eucDist(c2, c1) != sqrt(8)) {
		return 2;
	}
	if (eucDist(c3, c3) != 0) {
		return 2;
	}
	return 0;
}

int testIsRight() {

	Coords c1(0.0, 0.0);
	Coords c2(2.0, 2.0);
	Coords c3(1.0, 0.0);
	Coords c4(-1, -1);

	if (isRight(c1, c2, c3) != 1) {
		std::cout << "isRight(c1, c2, c3)-->" << isRight(c1, c2, c3)
				<< std::endl;
		return 1;
	}
	if (isRight(c2, c1, c3) != -1) {
		return 2;
	}
	if (isRight(c1, c2, c2) != 0) {
		return 3;
	}
	if (isRight(c1, c2, c1) != 0) {
		return 4;
	}
	if (isRight(c2, c2, c2) != 0) {
		return 5;
	}
	if (isRight(c1, c2, c4) != 0) {
		return 6;
	}
	if (isRight(c4, c1, c2) != 0) {
		return 7;
	}
	return 0;
}




#endif /* TESTS_BASICTEST_H_ */
