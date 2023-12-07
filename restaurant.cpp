#include "main.h"

int MAXSIZE = 10;

struct customer
{
	/* private:
		int result;
		int ID;

	public:
		void setResult(int X) { result = X; }
		int getResult() { return result; }
		void setID(int ID) { this->ID = ID; }
		int getID() { return ID; } */
	int result;
	int ID;
};

// LAPSE
char caesar(char c, int step)
{
	if (!((c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z')))
		return '\0';
	char res = c + step;
	if (((step > ('Z' - c)) && (c >= 'A' && c <= 'Z')) || ((step > ('z' - c)) && (c >= 'a' && c <= 'z')))
	{
		if (c >= 'A' && c <= 'Z')
			res = 'A' + step - ('Z' - c) - 1;
		else
			res = 'a' + step - ('z' - c) - 1;
	}
	return res;
}

string encrypt(string s)
{
	if (s.length() == 0)
		return s;
	unordered_map<char, int> M;
	for (int i = 0; i < s.length(); i++)
		M[s[i]]++;
	vector<pair<char, int>> freq(M.begin(), M.end());
	unordered_map<char, int> newM;
	for (auto &i : freq)
	{
		char newChar = caesar(i.first, i.second);
		newM[newChar] += i.second;
	}

	freq = vector<pair<char, int>>(newM.begin(), newM.end());
	sort(
		freq.begin(), freq.end(), [&](const pair<char, int> &a, const pair<char, int> &b)
		{
            if (a.second < b.second)
                return true;
            else if (a.second == b.second && (s.find(a.first) < s.find(b.first)))
                return true;
            else
                return false; });
	string result = "";
	for (auto i : freq)
		result += i.first;
	return result;
}

int decrypt(string s)
{
	if (s.length() == 0)
		return 0;
	if (s.length() > 10)
		s = s.substr(s.length() - 10);
	int result = 0;
	for (int i = s.length() - 1; i >= 0; i--)
	{
		if (s[i] == '1')
			result += pow(2, s.length() - 1 - i);
	}
	return result;
}
// Huffman : do later

enum res
{
	G = 0,
	S
};

// BST
template <class T>
class BST
{
protected:
	struct Node
	{
		T data;
		Node *left;
		Node *right;
		Node() : left(nullptr), right(nullptr) {}
		Node(const T &data, Node *left = nullptr, Node *right = nullptr) : data(data), left(left), right(right) {}
		Node(T &&data, Node *left = nullptr, Node *right = nullptr) : data(data), left(left), right(right) {}
	};
	Node *root;
	vector<T> postOrder;
	void clear(Node *root)
	{
		if (root == nullptr)
			return;
		else
		{
			clear(root->left);
			clear(root->right);
			delete root;
		}
	}
	void insert(const T &data, Node *&root)
	{
		if (root == nullptr)
			root = new Node(data);
		else
		{
			if (data < root->data)
				insert(data, root->left);
			else
				insert(data, root->right);
		}
	}
	T *find(const T &data, Node *root)
	{
		if (root)
			return root->data == data ? &root->data : (root->data > data ? find(data, root->left) : find(data, root->right));
		else
			return nullptr;
	}
	void remove(Node *&root, T value)
	{
		if (root == nullptr)
			return;
		if (value < root->value)
			remove(root->left, value);
		else if (value > root->value)
			remove(root->right, value);
		else
		{
			if (!root->left && !root->right)
			{
				delete root;
				root = nullptr;
			}
			else if (root->left && !root->right)
			{
				Node *temp = root;
				root = root->left;
				delete temp;
			}
			else if (root->right && !root->left)
			{
				Node *temp = root;
				root = root->right;
				delete temp;
			}
			else
			{
				Node *minRight = root->right;
				while (minRight->left != nullptr)
					minRight = minRight->left;
				root->value = minRight->value;
				remove(root->right, minRight->value);
			}
		}
	}
	void setPostOrder(Node *root)
	{
		if (root)
		{
			postOrder(root->left);
			postOrder(root->right);
			postOrder.push_back(root->data);
		}
	}

public:
	BST() : root(nullptr) {}
	~BST() { this->clear(); }
	void clear()
	{
		this->clear(this->root);
	}

	void insert(const T &data)
	{
		// insert(data, this->root);
		Node **p = &root;
		while (*p)
		{
			if ((*p)->data > data)
				p = &((*p)->left);
			else
				p = &((*p)->right);
		}
		*p = new Node(data);
	}
	void remove(const T &data)
	{
		remove(this->root, data);
	}
	T *find(const T &data)
	{
		// return find(data,root);
		Node *p = this->root;
		while (p)
		{
			if (p->data == data)
				return &p->data;
			else
				p = p->data > data ? p->left : p->right;
		}
		return nullptr;
	}
	vector<T> getPostOrder() { return postOrder; }
};

int numOfWays(vector<int> postOrder)
{
	// do later
	return 0;
}
// maybe I don't need this
/* class resG
{
private:
	vector<areaG> area;

public:
	void add(int ID, int val)
	{
		area[ID - 1].addCus(val);
	}
	void remove()
	{
		// For example the array (vector) is area, then call area[i]->removeCus() for all i
		for (int i = 0; i < MAXSIZE; i++)
		{
			area[i].removeCus();
		}
	}
}; */

class areaG
{
private:
	queue<customer> cusOrder;
	BST<int> *root;

public:
	areaG() : root(nullptr) {}
	void addCus(customer cus)
	{
		root->insert(cus.result);
		cusOrder.push(cus);
	}
	void removeCus()
	{
		if (root == nullptr)
			return;
		int k = numOfWays(root->getPostOrder());
		if (k >= cusOrder.size())
		{
			root->clear();
			while (!cusOrder.empty())
				cusOrder.pop();
		}
		else
		{
			for (int i = 0; i < k; i++)
			{
				root->remove(cusOrder.front().result);
				cusOrder.pop();
			}
		}
	}
	~areaG()
	{
		root->clear();
		delete root;
	}
};

class areaS
{
private:
	queue<customer> cusOrder;
	int size; // value
	int ID;

public:
	areaS(int ID) : size(0), ID(ID) {}
	void setID(int ID) { this->ID = ID; }
	int getID() { return ID; }
	void setSize(int size) { this->size = size; }
	int getSize() { return size; }
	void addCus(customer cus)
	{
		size++;
		cusOrder.push(cus);
	}
	void removeCus(int amount)
	{
		while (!cusOrder.empty() || amount > 0)
		{
			cusOrder.pop();
			size--;
			amount--;
		}
	}
	// a print method
	~areaS();
};

class resS
{ // a min heap
private:
	vector<areaS> area;
	vector<int> areaOrder; // use as queue to know which area comes first
	// store most recent updated
	vector<int> recentUsed; // the first element is the most recent, last element is the least recent

public:
	resS();

	int findRecentArea(int ID)
	{
		int idx = 0;
		while (recentUsed[idx] != ID)
			idx++;
		return idx;
	}

	int findOrderArea(int ID)
	{
		int idx = 0;
		while (areaOrder[idx] != ID)
			idx++;
		return idx;
	}

	int findArea(int ID)
	{
		int idx = 0;
		while (area[idx].getID() != ID)
			idx++;
		return idx;
	}

	void updateRecent(int ID)
	{
		// find area ID index
		int idx = findRecentArea(ID);
		// just leave it here, will modify later for the right purpose
		// in case that area has no customers left
		/* if (area[idx].getSize() == 0)
		{
			recentUsed.erase(recentUsed.begin() + idx);
			return;
		} */
		while (idx > 0)
		{
			swap(recentUsed[idx], recentUsed[idx - 1]);
			idx--;
		}
	}

	void updateOrder(int ID)
	{
		//
		areaOrder.push_back(ID);
	}

	void swapArea(areaS &a, areaS &b)
	{
		int tempSize = a.getSize();
		int tempID = a.getID();
		a.setSize(b.getSize());
		a.setID(b.getID());
		b.setSize(tempSize);
		b.setID(tempID);
	}
	// heap starts with 0
	/* void heapUp(int idx)
	{
		while (idx > 0)
		{
			int pIdx = idx / 2; // parent
			if (area[idx].getSize() >= area[pIdx].getSize() || (area[idx].getSize() == area[pIdx].getSize() && findOrderArea(area[idx].getID()) >= findOrderArea(area[pIdx].getID())))
				return;
			// swap
			swapArea(area[idx], area[pIdx]);
			int tempSize = area[idx].getSize();
			int tempID = area[idx].getID();
			area[idx].setSize(area[pIdx].getSize());
			area[idx].setID(area[pIdx].getID());
			area[pIdx].setSize(tempSize);
			area[pIdx].setID(tempID);
			idx = pIdx;
		}
	} */

	void heapDown(int idx)
	{
		//
		while (idx < area.size())
		{
			int cIdx = idx * 2 + 1; // left child
			if (cIdx >= area.size())
				return;
			if (cIdx + 1 < area.size()) // see which child will be used to swap
			{
				if (area[cIdx + 1].getSize() > area[cIdx].getSize() || (area[cIdx + 1].getSize() == area[cIdx].getSize() && findOrderArea(area[cIdx + 1].getID()) < findOrderArea(area[cIdx].getID())))
					cIdx++; // take right child if it is greater
			}
			if (area[idx].getSize() <= area[cIdx].getSize() || (area[idx].getSize() == area[cIdx].getSize() && findOrderArea(area[idx].getID()) <= findOrderArea(area[cIdx].getID())))
				return;
			else
			{
				swapArea(area[idx], area[cIdx]);
				idx = cIdx;
			}
		}
	}

	void buildHeap()
	{
		for (int i = area.size() / 2; i > -1; --i)
			heapDown(i);
	}

	int removeArea(int ID)
	{
		// set last elements to the removing element
		int idx = findArea(ID);
		area[idx].setSize(area[area.size() - 1].getSize());
		area[idx].setID(area[area.size() - 1].getID());
		area.pop_back();
		// do heapDown here
		heapDown(idx);
		// remove in areaOrder
		areaOrder.erase(areaOrder.begin() + findOrderArea(ID));
		// remove in recentUsed
		recentUsed.erase(recentUsed.begin() + findRecentArea(ID));
	}

	void remove(int NUM) // keiteiken
	{					 // heap down

		// choose NUM areas to remove NUM customers
		// find NUM area
		vector<int> removalArea; // contains areaS ID
		int idx = 0;
		while (removalArea.size() < NUM)
		{
			// TODO
		}
		// got all ID for deleting customers
		for (auto i : removalArea)
		{
			// update recentUsed
			updateRecent(i);
			// perform deleting
			int idx = findArea(i);
			area[idx].removeCus(NUM);
		}
		// after deleting in NUM areas
		// find if any area has 0 size
		for (auto i : area)
		{
			if (i.getSize() == 0)
				removeArea(i.getID());
		}
		// buildHeap again
		buildHeap();
	}

	void addArea(areaS newArea)
	{
		area.push_back(newArea);
		updateRecent(newArea.getID());
		updateOrder(newArea.getID());
		buildHeap();
	}

	void addCus(int ID, customer cus)
	{
		int index = findArea(ID);
		area[index].addCus(cus);
		updateRecent(ID);
		buildHeap();
	}
	~resS();
};

// Ingame function
void LAPSE(string name); // ecrypt and choose restaurant for customer
void KOKUSEN();			 // remove customers in resG
void KEITEIKEN(int NUM); // remove customers in resS
void HAND();			 // print Huffman tree
void LIMITLESS(int NUM); // print BST at area NUM in resG
void CLEAVE(int NUM);	 // print NUM customers info (LIFO)??? of areas by pre-order.

void simulate(string filename)
{
	vector<areaG> resG(MAXSIZE); // fixed size
	// vector<areaS> resS;			 // unfixed size <= MAXSIZE

	cout << "Good Luck";
	return;
}