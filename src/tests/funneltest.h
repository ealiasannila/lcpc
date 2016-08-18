/*
 * funneltest.h
 *
 *  Created on: Aug 17, 2016
 *      Author: elias
 */

#ifndef TESTS_FUNNELTEST_H_
#define TESTS_FUNNELTEST_H_
#include <vector>

void printChain(std::deque<Coords*> chain) {
	std::cout << "chain size: " << chain.size() << std::endl;
	for (Coords* c : chain) {
		std::cout << c->toString() << " <- " << c << std::endl;
	}
	std::cout << "\n";
}

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
	fun.reactToOpposite(&cl2);
	if (fun.getBase().getL() != &cl2 or fun.getBase().getR() != &cr1) {
		return 1;
	}
	fun.reactToOpposite(&cr2);
	if (fun.getBase().getL() != &cl2 or fun.getBase().getR() != &cr2) {
		return 2;
	}
	return 0;

}

int funnelGetChannelsTest() {
	Coords ca(0.0, 0.0);
	Coords cl1(0.5, 1.0);
	Coords cl2(0.4, 2.0);
	Coords cr1(1.0, 0.5);
	Coords cr2(2.0, 0.4);
	Coords s(3, 3);

	Funnel fun(&cl1, &ca, &cr1);
	fun.reactToOpposite(&cl2);
	fun.reactToOpposite(&cr2);
	std::deque<Coords*> funlc = fun.getLC();
	std::deque<Coords*> funrc = fun.getRC();
	if (funlc[0] != & ca or funlc[1] != &cl1 or funlc[2] != &cl2) {
		printChain(funlc);
		std::cout << "Left channel is not right" << std::endl;
		return 1;
	}

	if (funrc[0] != &ca or funrc[1] != &cr1 or funrc[2] != &cr2) {
		std::cout << "Right channel is not right" << std::endl;
		return 1;
	}
	return 0;
}

int funnelSplitTest() {
	Coords ca(0.0, 0.0);
	Coords cl1(0.5, 1.0);
	Coords cl2(0.4, 2.0);
	Coords cr1(1.0, 0.5);
	Coords cr2(2.0, 0.4);
	Coords s(3, 3);

	Funnel fun(&cl1, &ca, &cr1);
	fun.reactToOpposite(&cl2);
	fun.reactToOpposite(&cr2);
	Funnel f = fun.split(&s);
	/*std::cout << "--FUN:--" << std::endl;
	 std::cout << fun.toString() << std::endl;
	 std::cout << "--F:--\n" << std::endl;
	 std::cout << f.toString() << std::endl;
	 */
	std::deque<Coords*> funlc = fun.getLC();
	if (funlc[0] != &ca or funlc[1] != &s or funlc.size() != 2) {
		std::cout << "Left channel is not right after split in original funnel"
				<< std::endl;
		return 1;
	}
	std::deque<Coords*> funrc = fun.getRC();
	if (funrc[0] != &ca or funrc[1] != &cr1 or funrc[2] != &cr2 or funrc.size() != 3) {
		std::cout << "Right channel is not right after split in original funnel"
				<< std::endl;
		return 1;
	}

	std::deque<Coords*> frc = f.getRC();
	if (frc[0] != &ca or frc[1] != &s or frc.size() != 2) {
		std::cout << "Right channel is not right after split in new funnel"
				<< std::endl;
		return 1;
	}
	std::deque<Coords*> flc = f.getLC();
	if (flc[0] != &ca or flc[1] != &cl1 or flc[2] != &cl2 or flc.size() != 3) {
		std::cout << "Left channel is not right after split in new funnel"
				<< std::endl;
		return 1;
	}

	return 0;
}

int funnelShrinkTest() {
	Coords ca(0.0, 0.0);
	Coords cl1(0.5, 1);
	Coords cl2(0.7, 2);
	Coords cl3(0.7, 3);
	Coords ls(0.71, 2);

	Coords cr1(1.0, 0.5);
	Coords cr2(2.0, 0.7);
	Coords cr3(3.0, 0.7);
	Coords rs(2, 0.71);

	Funnel fun(&cl1, &ca, &cr1);

	fun.reactToOpposite(&cr2);
	fun.reactToOpposite(&cr3);
	fun.reactToOpposite(&cl2);
	fun.reactToOpposite(&cl3);
	fun.reactToOpposite(&rs);
	std::deque<Coords*> rc = fun.getRC();
	std::deque<Coords*> lc = fun.getLC();

	if (rc[0] != &ca or rc[1] != &cr1 or rc[2] != &rs or rc.size() != 3) {
		printChain(rc);
		std::cout << "Right Chain not correct after shrinkR" << std::endl;
		return 1;
	}
	if (lc[0] != &ca or lc[1] != &cl1 or lc[2] != &cl2 or lc[3] != &cl3 or lc.size() != 4) {
		std::cout << "Left chain not correct after shirnkR" << std::endl;
		return 2;
	}
	fun.reactToOpposite(&ls);
	lc = fun.getLC();
	if (lc[0] != &ca or lc[1] != &cl1 or  lc[2] != &ls or lc.size() != 3) {
		printChain(lc);
		std::cout << "Left chain not correct after shrinkL & shirnkR"
				<< std::endl;
		return 3;
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
	if (funnelGetChannelsTest() != 0) {
		return 1;
	}
	if (funnelSplitTest() != 0) {
		return 1;
	}
	if (funnelShrinkTest() != 0) {
		return 1;
	}
	return 0;
}

#endif /* TESTS_FUNNELTEST_H_ */
