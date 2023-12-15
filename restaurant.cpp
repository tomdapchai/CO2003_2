#include "main.h"

int MAXSIZE;

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

vector<pair<char, int>> encrypt(string &s)
{
	if (s.length() == 0)
	{
		vector<pair<char, int>> res;
		return res;
	}

	unordered_map<char, int> M;
	for (int i = 0; i < s.length(); i++)
		M[s[i]]++;
	vector<pair<char, int>> freq(M.begin(), M.end());
	unordered_map<char, int> newM;
	// combine same character
	for (auto &i : freq)
	{
		char newChar = caesar(i.first, i.second % 26);
		newM[newChar] += i.second;
	}
	// also encode s
	for (int i = 0; i < s.length(); i++)
	{
		s[i] = caesar(s[i], M[s[i]]);
	}
	freq = vector<pair<char, int>>(newM.begin(), newM.end());
	sort(
		freq.begin(), freq.end(), [&](const pair<char, int> &a, const pair<char, int> &b)
		{
            if (a.second < b.second)
                return true;
            else if (a.second == b.second) {
                if ((((a.first >= 'A' && a.first <= 'Z') && (b.first >= 'A' && b.first <= 'Z')) || ((a.first >= 'a' && a.first <= 'z') && (b.first >= 'a' && b.first <= 'z'))) && (a.first < b.first))
                    return true;
                else if ((a.first >= 'a' && a.first <= 'z') && (b.first >= 'A' && b.first <= 'Z'))
                    return true;
                else
                    return false;
            }
            else
                return false; });
	return freq;
}

int decrypt(string s)
{
	if (s.length() == 0)
		return 0;
	if (s.length() > 10)
		s = s.substr(s.length() - 10);
	int result = 0;
	for (int i = 0; i <= s.length() - 1; i++)
	{
		if (s[i] == '1')
			result += pow(2, i);
	}
	return result;
}
// Huffman : do later
struct HuffNode
{
	char ch;
	int freq;
	int order;
	HuffNode *left, *right;
};

// Function to allocate a new tree HuffNode
HuffNode *getHuffNode(char ch, int freq, HuffNode *left, HuffNode *right, int order)
{
	HuffNode *node = new HuffNode();

	node->ch = ch;
	node->freq = freq;
	node->left = left;
	node->right = right;
	node->order = order;
	return node;
}

// Comparison object to be used to order the heap
struct comp
{
	bool operator()(HuffNode *l, HuffNode *r)
	{
		// highest priority item has lowest frequency
		if (l->freq == r->freq)
			return l->order > r->order;
		return l->freq > r->freq;
	}
};

int height(HuffNode *root)
{
	int right = 0;
	int left = 0;
	if (root == NULL)
	{
		return 0;
	}
	// traverse until NULL -> temp = step, compare with result, swap if bigger, repeat until last node
	left = 1 + height(root->left);
	right = 1 + height(root->right);
	return left > right ? left : right;
}

// Function to get the balance factor of the node
int getBalanceFactor(HuffNode *node)
{
	if (node == nullptr)
		return 0;
	return height(node->left) - height(node->right);
}

// Function to perform a right rotation
HuffNode *rightRotate(HuffNode *root)
{
	HuffNode *temp = root->left;
	root->left = temp->right;
	temp->right = root;
	return temp;
}

// Function to perform a left rotation
HuffNode *leftRotate(HuffNode *root)
{
	HuffNode *temp = root->right;
	root->right = temp->left;
	temp->left = root;
	return temp;
}

// Builds Huffman Tree and decode given input text
bool isBalanced(HuffNode *root)
{
	return abs(getBalanceFactor(root)) <= 1;
}

HuffNode *balanceTree(HuffNode *root, int &counter)
{
	if (root == nullptr || isBalanced(root) || counter >= 3)
		return root;
	int balanceFactor = getBalanceFactor(root);
	if (balanceFactor > 1)
	{
		if (getBalanceFactor(root->left) >= 0)
		{
			root = rightRotate(root);
		}
		else
		{
			root->left = leftRotate(root->left);
			root = rightRotate(root);
		}
	}
	if (balanceFactor < -1)
	{
		if (getBalanceFactor(root->right) <= 0)
		{
			root = leftRotate(root);
		}
		else
		{
			root->right = rightRotate(root->right);
			root = leftRotate(root);
		}
	}
	counter++;
	// root->left = balanceTree(root->left, counter);
	// root->right = balanceTree(root->right, counter);
	return root;
}

// Function to traverse the tree in PreOrder
void preOrder(HuffNode *&root, HuffNode *&initialRoot, int &counter)
{
	if (root == nullptr || counter >= 3)
		return;

	bool rotated = false;
	if (!isBalanced(root))
	{
		root = balanceTree(root, counter);
		rotated = true;
	}

	// rotate intialRoot again if rotated is true
	if (rotated)
	{
		while (!isBalanced(initialRoot) && counter < 3)
			initialRoot = balanceTree(initialRoot, counter);
		preOrder(initialRoot, initialRoot, counter);
	}
	else
	{
		preOrder(root->left, initialRoot, counter);
		preOrder(root->right, initialRoot, counter);
	}
}

// Function to balance the entire tree
HuffNode *balanceEntireTree(HuffNode *root)
{
	int counter = 0;
	while (!isBalanced(root) && counter < 3)
		root = balanceTree(root, counter);

	if (isBalanced(root) && counter < 3)
		preOrder(root, root, counter);
	return root;
}

// Builds Huffman Tree and decode given input text
HuffNode *buildHuffmanTree(vector<pair<char, int>> freq)
{
	priority_queue<HuffNode *, vector<HuffNode *>, comp> pq;

	// Create a leaf HuffNode for each character and add it
	// to the priority queue.
	int order = 0;
	for (auto pair : freq)
	{
		// cout << pair.first << ": " << pair.second << endl;
		pq.push(getHuffNode(pair.first, pair.second, nullptr, nullptr, order));
		order++;
	}

	// do till there is more than one HuffNode in the queue
	while (pq.size() != 1)
	{
		// Remove the two HuffNodes of highest priority
		// (lowest frequency) from the queue
		HuffNode *left = pq.top();
		pq.pop();
		HuffNode *right = pq.top();
		pq.pop();

		// Create a new internal HuffNode with these two HuffNodes
		// as children and with frequency equal to the sum
		// of the two HuffNodes' frequencies. Add the new HuffNode
		// to the priority queue.
		int sum = left->freq + right->freq;
		HuffNode *node = getHuffNode('\0', sum, left, right, order);
		order++;

		// Balance the tree
		// node = balanceEntireTree(node);
		/* int counter = 0;
		node = balanceTree(node, counter); */
		node = balanceEntireTree(node);

		pq.push(node);
	}

	// root stores pointer to root of Huffman Tree
	HuffNode *root = pq.top();
	return root;
}

void printHuffTree(HuffNode *root)
{
	if (root == nullptr)
		return;

	printHuffTree(root->left);
	if (root->ch == '\0')
		cout << root->freq << "\n";
	else
		cout << root->ch << "\n";
	printHuffTree(root->right);
}

// Function to encode the Huffman tree
void encode(HuffNode *root, string str,
			unordered_map<char, string> &huffmanCode)
{
	if (root == nullptr)
		return;

	// found a leaf node
	if (!root->left && !root->right)
	{
		huffmanCode[root->ch] = str;
	}

	encode(root->left, str + "0", huffmanCode);
	encode(root->right, str + "1", huffmanCode);
}

// Function to build and print Huffman Codes
string getCode(HuffNode *root, string name)
{
	unordered_map<char, string> huffmanCode;
	encode(root, "", huffmanCode);

	string result = "";
	/* for (auto i : huffmanCode)
	{
		cout << i.first << ": " << i.second << endl;
		// result = result + i.second + "\n";
	} */
	for (int i = 0; i < name.length(); i++)
	{
		result += huffmanCode[name[i]];
	}
	return result;
}

// Function to free the Huffman tree
void freeHuffman(HuffNode *&root)
{
	if (root == nullptr)
		return;
	freeHuffman(root->left);
	freeHuffman(root->right);
	delete root;
	root = nullptr;
}

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
		if (value < root->data)
			remove(root->left, value);
		else if (value > root->data)
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
				root->data = minRight->data;
				remove(root->right, minRight->data);
			}
		}
	}
	void setPostOrder(Node *root)
	{
		if (root)
		{
			setPostOrder(root->left);
			setPostOrder(root->right);
			postOrder.push_back(root->data);
		}
	}

	void inOrder(Node *root)
	{
		if (root)
		{
			inOrder(root->left);
			cout << root->data << endl;
			inOrder(root->right);
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
	vector<T> getPostOrder()
	{
		setPostOrder(root);
		vector<int> copy(postOrder.begin(), postOrder.end());
		postOrder.clear();
		return copy;
	}
	void printInOrder()
	{
		inOrder(root);
	}
};

// modified version of factorial
long long int factorial(int n, int start = 1)
{
	long long int fact = 1;
	for (int i = start; i <= n; i++)
		fact = (fact * i);
	return fact;
}

// calculate number of ways
long long int numOfWays(vector<int> postOrder)
{
	if (postOrder.size() <= 1)
		return 1;

	vector<int> LST, RST;
	int root = postOrder.back();
	//  Splitting into left and right subtrees
	for (int i = 0; i < postOrder.size() - 1; i++)
	{
		if (postOrder[i] < root)
			LST.push_back(postOrder[i]);
		else
			RST.push_back(postOrder[i]);
	}
	// Recursive calls for left and right subtrees
	long long int left = numOfWays(LST);
	long long int right = numOfWays(RST);
	// Calculating number of ways using formula
	int big = LST.size() > RST.size() ? LST.size() : RST.size(); // bigger sub tree
	int small = RST.size() + LST.size() - big;					 // smaller sub tree
	// find which tree has bigger size to perform reducing factorial:
	/* instead of (big + small)! / (big! * small!), split it into:
		((big + small)! / big!) / small!
		then it equals to:
		(big + 1) * (big + 2) *...* (big + small) / small!
	*/
	long long int ways = factorial(big + small, big + 1) / factorial(small);
	return ways * left * right;
}
class areaG
{
private:
	queue<int> cusOrder;
	BST<int> *tree;

public:
	areaG()
	{
		tree = new BST<int>;
	}
	void addCus(int cus)
	{
		tree->insert(cus);
		cusOrder.push(cus);
		// cout << cusOrder.size() << endl;
	}
	void removeCus()
	{
		if (tree == nullptr || cusOrder.empty())
			return;
		vector<int> pO = tree->getPostOrder();
		/* for (int i : pO)
			cout << i << " "; */
		int k = numOfWays(pO) % MAXSIZE;
		// cout << k << endl;
		//  cout << "size " << cusOrder.size() << endl;
		if (k > 0)
		{
			if (k >= cusOrder.size())
			{
				tree->clear();
				tree = new BST<int>;
				while (!cusOrder.empty())
					cusOrder.pop();
				// cout << "true\n";
			}
			else
			{
				for (int i = 0; i < k; i++)
				{
					tree->remove(cusOrder.front());
					cusOrder.pop();
				}
			}
		}
	}
	void print()
	{
		// cout << "size: " << cusOrder.size() << endl;
		if (cusOrder.empty())
		{
			cout << "empty\n";
			return;
		}
		tree->printInOrder();
	}
	// BST<int> *getTree() { return tree; }
	~areaG()
	{
		tree->clear();
		delete tree;
	}
};

class areaS
{
private:
	queue<int> cusOrder;
	int size; // value
	int ID;

public:
	areaS(int ID) : size(0), ID(ID) {}
	void setID(int ID) { this->ID = ID; }
	int getID() { return ID; }
	void setSize(int size) { this->size = size; }
	int getSize() { return size; }
	void setOrder(queue<int> order) { cusOrder = order; }
	queue<int> getOrder() { return cusOrder; }
	void addCus(int cus)
	{
		size++;
		cusOrder.push(cus);
	}
	void removeCus(int amount)
	{
		while (!cusOrder.empty() && amount > 0)
		{
			cout << cusOrder.front() << "-" << ID << endl;
			cusOrder.pop();
			size--;
			amount--;
		}
	}
	void printInfo(int NUM)
	{
		if (cusOrder.empty())
			return;
		stack<int> temp1;
		stack<int> temp2;
		// throw to stack
		while (!cusOrder.empty())
		{
			temp1.push(cusOrder.front());
			cusOrder.pop();
		}
		// print first NUM values in stack
		while (!temp1.empty())
		{
			if (NUM > 0)
				cout << ID << "-" << temp1.top() << endl;
			temp2.push(temp1.top());
			temp1.pop();
			NUM--;
		}
		// restore queue
		while (!temp2.empty())
		{
			cusOrder.push(temp2.top());
			temp2.pop();
		}
	}
	~areaS() {}
};

void swapArea(areaS &a, areaS &b)
{
	int tempSize = a.getSize();
	int tempID = a.getID();
	queue<int> tempOrder = a.getOrder();
	a.setOrder(b.getOrder());
	a.setSize(b.getSize());
	a.setID(b.getID());
	b.setSize(tempSize);
	b.setID(tempID);
	b.setOrder(tempOrder);
}

class resS
{ // a min heap
private:
	vector<areaS> area;
	vector<int> recentUsed;
	// vector<int> areaOrder; // use as queue to know which area comes first
	//  store most recent updated
	// the first element is the most recent, last element is the least recent
	void addArea(areaS newArea)
	{
		area.push_back(newArea);
		updateRecent(newArea.getID());
		buildHeap();
	}

	void removeArea(int ID)
	{
		// set last elements to the removing element
		int idx = findArea(ID);
		area[idx].setSize(area[area.size() - 1].getSize());
		area[idx].setID(area[area.size() - 1].getID());
		area[idx].setOrder(area[area.size() - 1].getOrder());
		area.pop_back();
		// do heapDown here
		heapDown(idx);
		// remove in recentUsed
		recentUsed.erase(recentUsed.begin() + findRecentArea(ID));
	}

public:
	resS() {}

	int findRecentArea(int ID)
	{
		if (recentUsed.size() == 0)
			return -1;
		int idx = 0;
		while (idx < recentUsed.size() && recentUsed[idx] != ID)
			idx++;
		return idx < recentUsed.size() ? idx : -1;
	}

	int findArea(int ID)
	{
		int idx = 0;
		while (idx < area.size() && area[idx].getID() != ID)
			idx++;
		return idx < area.size() ? idx : -1;
	}

	void updateRecent(int ID) // the most recent used area is at the begin, least is at the end
	{
		// find area ID index
		int idx = findRecentArea(ID);
		// If this is newly added area
		if (idx == -1)
		{
			recentUsed.push_back(ID);
			idx = recentUsed.size() - 1;
		}
		while (idx > 0)
		{
			swap(recentUsed[idx], recentUsed[idx - 1]);
			idx--;
		}
	}
	// heap starts with 0
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
				if (area[cIdx + 1].getSize() < area[cIdx].getSize() || (area[cIdx + 1].getSize() == area[cIdx].getSize() && findRecentArea(area[cIdx + 1].getID()) > findRecentArea(area[cIdx].getID())))
					cIdx++; // take right child if it is smaller (if equal size choose which is less recently used)
			}
			if (area[idx].getSize() < area[cIdx].getSize() || (area[idx].getSize() == area[cIdx].getSize() && findRecentArea(area[idx].getID()) > findRecentArea(area[cIdx].getID())))
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

	void printPreOrder(int NUM, int idx = 0)
	{
		if (idx < area.size())
		{
			area[idx].printInfo(NUM);
			printPreOrder(NUM, idx * 2 + 1); // left child
			printPreOrder(NUM, idx * 2 + 2); // right child
		}
	}

	void remove(int NUM) // keiteiken
	{
		// choose NUM areas to remove NUM customers
		// find NUM area
		vector<areaS> copyArea(area.begin(), area.end());
		// sort base on size and least used
		sort(copyArea.begin(), copyArea.end(), [&](areaS &a, areaS &b)
			 {
			if (a.getSize() < b.getSize())
				return true;
			else if (a.getSize() == b.getSize()) {
				// If a is less recent used than b then true
				if (findRecentArea(a.getID()) > findRecentArea(b.getID()))
					return true;
				else
					return false;
			}
			else
				return false; });
		vector<int> removalArea; // contains ID of NUM areaS to be removed

		// Load first NUM areas in copyArea to removalArea, if NUM > area.size() then print all area
		for (int i = 0; i < NUM && i < copyArea.size(); i++)
		{
			removalArea.push_back(copyArea[i].getID());
		}
		// got all area ID for deleting customers
		for (auto &i : removalArea)
		{
			// perform deleting
			int idx = findArea(i);
			area[idx].removeCus(NUM);
			// after deleting if no customer left, remove area
			if (area[idx].getSize() <= 0)
				removeArea(i);
			else
				updateRecent(i);
			// reheap
			buildHeap();
		}
	}

	void addCus(int ID, int cus)
	{
		int index = findArea(ID);
		if (index != -1)
		{
			area[index].addCus(cus);
			updateRecent(ID);
			buildHeap();
		}
		else // If the area hasn't existed, create new one
		{
			areaS newArea(ID);
			newArea.addCus(cus);
			addArea(newArea);
		}
	}
	~resS()
	{
		area.clear();
		recentUsed.clear();
	}
};

// Ingame function

void simulate(string filename)
{
	vector<areaG> ResG; // fixed size
	resS ResS;			// unfixed size <= MAXSIZE
	ifstream ss(filename);
	string str, maxsize, name, num;
	HuffNode *root = nullptr;
	while (ss >> str)
	{
		if (str == "MAXSIZE")
		{
			ss >> maxsize;
			MAXSIZE = stoi(maxsize);
			ResG.resize(MAXSIZE);
		}
		else if (str == "LAPSE")
		{
			ss >> name;
			if (name.length() < 3)
				continue;
			vector<pair<char, int>> freq = encrypt(name);
			if (freq[0].second == name.length())
				continue;
			freeHuffman(root);
			root = buildHuffmanTree(freq);
			int result = decrypt(getCode(root, name));
			int ID = result % MAXSIZE + 1;
			// cout << "ID: " << ID << endl;
			if (result % 2 == 1) // Gojo
			{
				// BST starts from 0
				ResG[ID - 1].addCus(result);
			}
			else // Sukuna
			{
				ResS.addCus(ID, result);
			}
			//
		}
		else if (str == "KOKUSEN")
		{
			for (areaG &area : ResG)
			{
				area.removeCus();
			}
		}
		else if (str == "KEITEIKEN")
		{
			ss >> num;
			int NUM = stoi(num);
			ResS.remove(NUM);
		}
		else if (str == "HAND")
		{
			// print Huffman
			// no idea what to print, but gonna use the root
			printHuffTree(root);
		}
		else if (str == "LIMITLESS")
		{
			// cout << "LIMITLESS\n";
			ss >> num;
			int NUM = stoi(num);
			cout << NUM << endl;
			if (NUM <= MAXSIZE && NUM >= 1)
				ResG[NUM - 1].print();
		}
		else if (str == "CLEAVE")
		{
			// cout << "CLEAVE\n";
			ss >> num;
			int NUM = stoi(num);
			ResS.printPreOrder(NUM);
		}
		else if (str == "STOP")
		{
			cout << "STOP\n";
			break;
		}
	}
	ss.close();
	// add a function to free the Huffman to avoid memleak
	freeHuffman(root);
	// cout << "Good Luck";
	return;
}