/*
 * UnionFind.h
 *
 *  Created on: Jan 29, 2012
 *      Author: aguthals
 */

#ifndef UNIONFIND_H_
#define UNIONFIND_H_

//TR1 includes. GCC 4.0 and above only!
#ifdef __GLIBCXX__
#  include <tr1/unordered_map>
#else
#  ifdef __IBMCPP__
#    define __IBMCPP_TR1__
#  endif
#  include <unordered_map>
#endif
#include <vector>

using namespace std;

namespace specnets {

template<typename T> struct UnionFindNode {
	UnionFindNode<T>* parent;
	T data;
};

template<typename T> class UnionFind {

protected:
	typedef tr1::unordered_map<T, size_t> NodeMap;
	typedef vector<UnionFindNode<T> > NodeContainer;

	NodeContainer* nodeContainer;
	NodeMap* nodeMap;

public:
	UnionFind() :
		nodeContainer(0x0), nodeMap(0x0) {
	}

	template<typename InputIterator> UnionFind(InputIterator first,
			InputIterator last) :
		nodeContainer(0x0), nodeMap(0x0) {
		initialize(first, last);
	}

	UnionFind(vector<T>& inputData) :
		nodeContainer(0x0), nodeMap(0x0) {
		initialize(inputData);
	}

	~UnionFind() {
		if (nodeContainer) {
			delete nodeContainer;
		}
		if (nodeMap) {
			delete nodeMap;
		}
	}

	template<typename InputIterator> void initialize(InputIterator first,
			InputIterator last) {
		if (nodeContainer) {
			nodeContainer->resize(0);
		} else {
			nodeContainer = new NodeContainer(0);
		}

		int curIdx = 0;
		UnionFindNode<T> tempNode;
		tempNode.parent = 0;
		while (first != last) {
			tempNode.data = *first;
			nodeContainer->push_back(tempNode);
			first++;
			++curIdx;
		}

		if (nodeMap) {
			nodeMap->clear();
			nodeMap->rehash(nodeContainer->size());
		} else {
			nodeMap = new NodeMap(nodeContainer->size());
		}

		for (unsigned int i = 0; i < nodeContainer->size(); i++) {
			(*nodeMap)[(*nodeContainer)[i].data] = i;
		}
	}

	void initialize(vector<T>& inputData) {
		if (nodeContainer) {
			nodeContainer->resize(inputData.size());
		} else {
			nodeContainer = new NodeContainer(inputData.size());
		}

		for (unsigned int i = 0; i < inputData.size(); i++) {
			T* data = &inputData[i];
			(*nodeContainer)[i].parent = (UnionFindNode<T>*) 0;
			(*nodeContainer)[i].data = *data;
		}

		if (nodeMap) {
			nodeMap->clear();
			nodeMap->rehash(nodeContainer->size());
		} else {
			nodeMap = new NodeMap(nodeContainer->size());
		}

		for (unsigned int i = 0; i < nodeContainer->size(); i++) {
			(*nodeMap)[(*nodeContainer)[i].data] = i;
		}
	}

	UnionFindNode<T>* findParent(UnionFindNode<T>* inNode) {
		UnionFindNode<T>* parentNode = inNode;
		while (parentNode->parent) {
			parentNode = parentNode->parent;
		}
		if (parentNode != inNode) {
			inNode->parent = parentNode;
		}
		return parentNode;
	}

	UnionFindNode<T>* getNode(T item) {
		if (nodeMap->count(item) == 0) {
			ERROR_MSG("Could not find item " << item);
			return 0;
		}
		return &(*nodeContainer)[(*nodeMap)[item]];
	}

	void getSets(vector<vector<T> >& outputSets) {

		if (size() == 0) {
			outputSets.resize(0);
			return;
		}

		map<UnionFindNode<T>*, list<UnionFindNode<T>*> >* mappedSets = new map<
				UnionFindNode<T>*, list<UnionFindNode<T>*> > ();

		for (unsigned int i = 0; i < nodeContainer->size(); i++) {
			UnionFindNode<T>* child = &(*nodeContainer)[i];
			UnionFindNode<T>* parent = findParent(child);
			if (mappedSets->count(parent) == 0) {
				list<UnionFindNode<T>*> tempList;
				tempList.push_back(child);
				(*mappedSets)[parent] = tempList;
			} else {
				(*mappedSets)[parent].push_back(child);
			}
		}

		outputSets.resize(mappedSets->size());
		unsigned int idxUse = 0;

		for (typename map<UnionFindNode<T>*, list<UnionFindNode<T>*> >::iterator
				mapIt = mappedSets->begin(); mapIt != mappedSets->end(); mapIt++) {
			outputSets[idxUse].resize(mapIt->second.size());
			unsigned int j = 0;
			for (typename list<UnionFindNode<T>*>::iterator cIt =
					mapIt->second.begin(); cIt != mapIt->second.end(); cIt++) {
				outputSets[idxUse][j++] = (*cIt)->data;
			}
			++idxUse;
		}

		delete mappedSets;
	}

	bool unionSets(T set1, T set2) {
		UnionFindNode<T>* parent1 = findParent(getNode(set1));
		UnionFindNode<T>* parent2 = findParent(getNode(set2));

		if (parent1 == parent2) {
			return false;
		}
		parent2->parent = parent1;
		return true;
	}

	bool inSameSet(T set1, T set2) {
		UnionFindNode<T>* parent1 = findParent(getNode(set1));
		UnionFindNode<T>* parent2 = findParent(getNode(set2));

		return parent1 == parent2;
	}

	unsigned int size() {
		return nodeContainer->size();
	}

};

}

#endif /* UNIONFIND_H_ */
