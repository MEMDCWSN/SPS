/*
 * TrieMap.h
 *
 *  Created on: May 2, 2012
 *      Author: aguthals
 */

#ifndef TRIEMAP_H_
#define TRIEMAP_H_

#include "Logger.h"

#include <set>
#include <list>
#include <vector>

using namespace std;

namespace specnets
{

  template<class T> class TrieMap
  {
  private:

    static const unsigned int NUM_CHARS = 26;

    static unsigned int charToIdx(char character)
    {
      int c = (int)(toupper(character)) - (int)('A');

      if (c < 0 || c >= NUM_CHARS)
      {
        ERROR_MSG("Unknown character \'" << character << "\'!");
        abort();
      }
      return (unsigned int)c;
    }

    struct node
    {
      char character; // character of the node
      unsigned int eow; // number of complete words that end here
      unsigned int prefixes; // number of words that have this prefix
      set<T> data; // data referenced by this word
      vector<node*> children; // references to all possible children
      node* parent;
    }*root;

    node* getNode(const string& word, const T* valuePtr = NULL) const
    {
      if (word.length() == 0)
      {
        return NULL;
      }

      node *t = root;
      unsigned int i = 0;
      unsigned int wordLen = word.length();
      char s;
      unsigned int c;
      while (i < wordLen)
      {
        s = word.at(i);
        c = charToIdx(s);
        if (t->children[c] == NULL)
        {
          return NULL;
        }
        else
        {
          t = t->children[c];
        }
        i++;
      }

      if (t->eow > 0 && valuePtr == NULL)
      {
        return t;
      }
      else if (t->eow > 0 && valuePtr != NULL)
      {
        if (t->data.count(*valuePtr) > 0)
        {
          return t;
        }
        return NULL;
      }
      else
      {
        return NULL;
      }
    }

  public:

    TrieMap()
    {
      root = new node();
      root->character = '\0';
      root->eow = 0;
      root->prefixes = 0;
      root->data.clear();
      root->parent = NULL;
      root->children.resize(NUM_CHARS);
      for (unsigned int i = 0; i < NUM_CHARS; i++)
      {
        root->children[i] = NULL;
      }
    }

    ~TrieMap()
    {
      clear();
      delete root;
    }

    void clear()
    {
      list<node*> stack;
      for (unsigned int i = 0; i < NUM_CHARS; i++)
      {
        if (root->children[i] != NULL)
        {
          stack.push_back(root->children[i]);
        }
      }

      while (stack.size() > 0)
      {
        node* next = stack.back();
        stack.pop_back();
        for (unsigned int i = 0; i < NUM_CHARS; i++)
        {
          if (next->children[i] != NULL)
          {
            stack.push_back(next->children[i]);
            next->children[i]->parent = NULL;
            next->children[i] = NULL;
          }
        }
        next->data.clear();
        next->children.resize(0);
        delete next;
      }
      root->prefixes = 0;
      for (unsigned int i = 0; i < NUM_CHARS; i++)
      {
        root->children[i] = NULL;
      }
    }

    TrieMap<T>& operator=(const TrieMap<T>& other)
    {
      clear();

      list<pair<node*, node*> > queue;
      pair<node*, node*> next(root, other.root);
      pair<node*, node*> nextChild;
      queue.push_back(next);

      while (queue.size() > 0)
      {
        next = queue.front();
        queue.pop_front();

        next.first->character = next.second->character;
        next.first->prefixes = next.second->prefixes;
        next.first->eow = next.second->eow;
        next.first->data = next.second->data;

        for (unsigned int i = 0; i < NUM_CHARS; i++)
        {
          if (next.second->children[i] == NULL)
          {
            next.first->children[i] = NULL;
          }
          else
          {
            node* t = new node();

            t->data.clear();
            t->parent = next.first;
            t->children.resize(NUM_CHARS);
            next.first->children[i] = t;

            nextChild.first = t;
            nextChild.second = next.second->children[i];
            queue.push_back(nextChild);
          }
        }
      }
      return *this;
    }

    void insertWord(const string& word, T* valuePtr = NULL)
    {
      if (word.length() == 0)
      {
        return;
      }

      node *t = root;
      unsigned int i = 0;
      unsigned int wordLen = word.length();
      char s;
      unsigned int c;
      while (i < wordLen)
      {
        s = word.at(i);
        c = charToIdx(s);
        t->prefixes++;
        if (t->children[c] == NULL)
        {
          node* n = new node();
          n->character = s;
          n->eow = 0;
          n->prefixes = 0;
          n->data.clear();
          n->children.resize(NUM_CHARS);
          for (unsigned int j = 0; j < NUM_CHARS; j++)
          {
            n->children[j] = NULL;
          }
          t->children[c] = n;
          n->parent = t;
          t = n;
        }
        else
        {
          t = t->children[c];
        }
        i++;
      }
      t->eow++;
      t->prefixes++;
      if (valuePtr != NULL)
      {
        t->data.insert(*valuePtr);
      }
    }

    unsigned int size() const
    {
      return root->prefixes;
    }

    unsigned int findWord(const string& word, list<T>* foundValues = NULL) const
    {
      if (foundValues != NULL)
      {
        foundValues->clear();
      }
      node* n = getNode(word);

      if (n == NULL || n->eow == 0)
      {
        return 0;
      }

      if (foundValues != NULL)
      {
        foundValues->assign(n->data.begin(), n->data.end());
      }
      return n->eow;
    }

    unsigned int deleteWord(const string& word, T* valuePtr = NULL)
    {
      node* n = getNode(word, valuePtr);
      if (n == NULL || n->eow == 0)
      {
        return 0;
      }
      unsigned int numRem = 0;
      if (valuePtr != NULL)
      {
        n->data.erase(*valuePtr);
        numRem = 1;
      }
      else if (n->data.size() > 0)
      {
        numRem = n->data.size();
      }
      else
      {
        numRem = 1;
      }

      n->eow -= numRem;

      while (n->parent != NULL)
      {
        n->prefixes -= numRem;
        node* del = n;
        n = del->parent;
        if (del->prefixes == 0)
        {
          n->children[charToIdx(del->character)] = NULL;
          delete del;
        }
      }
      root->prefixes -= numRem;
      return numRem;
    }

    unsigned int deleteWords(const string& word)
    {
      node* n = getNode(word);
      if (n == NULL || n->eow == 0)
      {
        return 0;
      }
      unsigned int numRem = n->eow;

      n->eow -= numRem;

      while (n->parent != NULL)
      {
        n->prefixes -= numRem;
        node* del = n;
        n = del->parent;
        if (del->prefixes == 0)
        {
          n->children[charToIdx(del->character)] = NULL;
          delete del;
        }
      }
      root->prefixes -= numRem;
      return numRem;
    }

    void findPrefixes(const string& word, list<T>& outputPrefixes) const
    {
      outputPrefixes.clear();
      node* n = getNode(word);
      if (n == NULL)
      {
        return;
      }

      while (n->parent != NULL)
      {
        if (n->eow > 0)
        {
          outputPrefixes.insert(outputPrefixes.end(),
                                n->data.begin(),
                                n->data.end());
        }
        n = n->parent;
      }
    }

    bool replaceValue(const string& word, const T& oldValuePtr, const T& newValuePtr)
    {
      node* n = getNode(word, &oldValuePtr);

      if (n == NULL)
      {
        return false;
      }

      n->data.erase(oldValuePtr);
      n->data.insert(newValuePtr);
      return true;
    }

  };

  typedef TrieMap<void*> Trie;
}

#endif /* TRIEMAP_H_ */
