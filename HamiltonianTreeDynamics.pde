int margin = 50;
int nodeRadius = 20;
double scaleY = 150;
int frame_rate = 50;

Tree tree;
double theta = 1.0;
int nLeaves = 5;

double dt = 0.02;

boolean paused = false;

void setup() {
  size(800, 600);
  frameRate(frame_rate);

  strokeWeight(4);
  stroke(#ffffff);
  fill(#0000ff);

  double[] leafAges = new double[nLeaves];
  for (int i=0; i<leafAges.length; i++)
    leafAges[i] = 0.1*i;
  tree = new Tree(leafAges, 1.0);
}

void draw() {
  if (!paused)
    tree.doStep();
  background(#000000);
  drawTree();
}

void keyPressed() {
  switch(key) {
    case ' ':
      paused = !paused;
      break;
    case 'r':
      setup();
      break;
    case '+':
      nLeaves += 1;
      setup();
      break;
    case '-':
      nLeaves = max(nLeaves-1, 3);
      setup();
      break;
    default:
      break;
  }
}

void mouseScrolled() {
  if (mouseScroll>0) {
    scaleY *= 1.1;
  } else {
    scaleY /= 1.1;
  }
}

void drawTree() {      

  for (Node node : tree.getNodeList ()) {
    if (node.isRoot())
      continue;

    int[] pos = node.getScreenPosition();
    int[] parentPos = node.getParent().getScreenPosition();
    line(pos[0], pos[1], parentPos[0], parentPos[1]);
  }

  for (Node node : tree.getNodeList ()) {
    int[] pos = node.getScreenPosition();
    ellipse(pos[0], pos[1], nodeRadius, nodeRadius);
  }
}

class Node {
  Node parent;
  ArrayList<Node> children;
  double age = 0;
  double layoutPosition=-1;
  int lineages;

  double momentum = 0.0;

  Node() {
    children = new ArrayList<Node>();
  }

  boolean isLeaf() {
    return children.isEmpty();
  }

  boolean isRoot() {
    return parent == null;
  }

  Node addChild(Node child) {
    children.add(child);
    child.parent = this;
    return this;
  }

  Node setParent(Node parent) {
    this.parent = parent;
    return this;
  }

  Node setAge(double age) {
    this.age = age;
    return this;
  }

  double getAge() {
    return age;
  }

  Node setMomentum(double momentum) {
    this.momentum = momentum;
    return this;
  }

  double getMomentum() {
    return momentum;
  }

  ArrayList<Node> getChildren() {
    return children;
  }

  Node getParent() {
    return parent;
  }

  int getLineages() {
    return lineages;
  }

  Node setLineages(int lineages) {
    this.lineages = lineages;
    return this;
  }

  ArrayList<Node> getChildList() {
    ArrayList<Node> res = new ArrayList<Node>();

    if (isLeaf()) {
      res.add(this);
    } else {
      for (int i=0; i<children.size (); i++) {
        res.addAll(children.get(i).getChildList());
      }
    }

    return res;
  }

  ArrayList<Node> getNodeList() {
    ArrayList<Node> res = new ArrayList<Node>();

    res.add(this);
    for (int i=0; i<children.size (); i++) {
      res.addAll(children.get(i).getNodeList());
    }

    return res;
  }

  double getLayoutPosition() {
    if (layoutPosition<0) {
      layoutPosition = 0.0;
      for (int i=0; i<children.size (); i++) {
        layoutPosition += children.get(i).getLayoutPosition();
      }
      layoutPosition /= children.size();
    }

    return layoutPosition;
  }

  String getNewick() {
    String newick = "";
    if (!children.isEmpty()) {
      newick += "(";
      for (int i=0; i<children.size (); i++) {
        if (i>0)
          newick += ",";
        newick += children.get(i).getNewick();
      }
      newick += ")";
    }

    newick += ":";
    if (!isRoot()) {
      newick += String.valueOf(parent.age - age);
    } else {
      newick += "0.0";
    }

    return newick;
  }

  int[] getScreenPosition() {
    int [] pos = new int[2];
    double scaleX = (width-2*margin)/(tree.getLeafCount()-1);
    pos[0] = margin + (int)(scaleX*getLayoutPosition());
    pos[1] = height - margin - (int)(scaleY*getAge());

    return pos;
  }
}

class Tree {
  Node root;
  ArrayList<Node> childList;
  ArrayList<Node> nodeList;

  double[] Q, P;

  // Create new tree with provided node as root.
  Tree(Node root) {
    this.root = root;
    childList = root.getChildList();
    nodeList = root.getNodeList();

    for (int i=0; i<childList.size (); i++) {
      childList.get(i).layoutPosition = i;
    }

    Q = new double[nodeList.size()];
    P = new double[nodeList.size()];
  }

  // Simulate a new coalescent tree with these leaf ages.
  Tree(double[] leafAges, double simTheta) {
    ArrayList<Node> inactiveNodes = new ArrayList<Node>();
    ArrayList<Node> activeNodes = new ArrayList<Node>();
    for (double age : leafAges) {
      inactiveNodes.add(new Node().setAge(age));
    }

    sortNodeList(inactiveNodes);

    double t = 0.0;
    while (activeNodes.size ()>1 || !inactiveNodes.isEmpty()) {

      int k = activeNodes.size();

      double propensity = 0.5*k*(k-1)/simTheta;

      boolean tIsInfinite = false;
      if (propensity > 0) {
        t += -log(random(1))/propensity;
      } else {
        tIsInfinite = true; // Hack to get around processing.js limitation
      }

      if (!inactiveNodes.isEmpty() && (t>inactiveNodes.get(0).getAge() || tIsInfinite)) {
        t = inactiveNodes.get(0).getAge();
        activeNodes.add(inactiveNodes.get(0));
        inactiveNodes.remove(0);
        continue;
      }

      // Select a pair of activeNodes to coalesce
      int idx1 = int(random(k));
      int idx2;
      do {
        idx2 = int(random(k));
      } 
      while (idx2 == idx1);

      Node node1 = activeNodes.get(idx1);
      Node node2 = activeNodes.get(idx2);

      Node newNode = new Node().setAge(t).addChild(node1).addChild(node2);
      activeNodes.remove(node1);
      activeNodes.remove(node2);
      activeNodes.add(newNode);
    }

    this.root = activeNodes.get(0);
    childList = root.getChildList();
    nodeList = root.getNodeList();

    for (int i=0; i<childList.size (); i++) {
      childList.get(i).layoutPosition = i;
    }

    Q = new double[nodeList.size()];
    P = new double[nodeList.size()];
  }

  void setRoot(Node root) {
    this.root = root;
  }

  ArrayList<Node> getChildList() {
    return childList;
  }

  ArrayList<Node> getNodeList() {
    return nodeList;
  }

  String getNewick() {
    return root.getNewick() + ";";
  }

  int getLeafCount() {
    return childList.size();
  }

  int getNodeCount() {
    return nodeList.size();
  }

  void sortNodeList(ArrayList<Node> nodeList) {
    boolean sorted;
    do {
      sorted = true;
      for (int i=0; i<nodeList.size ()-1; i++) {
        if (nodeList.get(i+1).getAge()<nodeList.get(i).getAge()) {
          Node tmp = nodeList.get(i+1);
          nodeList.set(i+1, nodeList.get(i));
          nodeList.set(i, tmp);
          sorted = false;
        }
      }
    } 
    while (!sorted);

    int k=0;
    for (Node node : nodeList) {
      if (node.isLeaf())
        k += 1;
      else
        k -= 1;
      node.setLineages(k);
    }
  }  

  void getQandP(double[] Q, double[] P) {
    sortNodeList(nodeList);
    for (int i=0; i<nodeList.size (); i++) {
      Q[i] = nodeList.get(i).getAge();
      P[i] = nodeList.get(i).getMomentum();
    }
  }

  void setQandP(double[] Q, double[] P) {
    for (int i=0; i<nodeList.size (); i++) {
      nodeList.get(i).setAge(Q[i]);
      nodeList.get(i).setMomentum(P[i]);
    }
  }

  void getQdotAndPdot(double[] Q, double[] P, double[] Qdot, double[] Pdot) {
    for (int i=0; i<nodeList.size (); i++) {
      if (nodeList.get(i).isLeaf()) {
        Qdot[i] = 0.0;
        Pdot[i] = 0.0;
        P[i] = 0.0;
        continue;
      } 

      int ki = nodeList.get(i).getLineages();
      int kiminus1 = nodeList.get(i-1).getLineages();
      Qdot[i] = P[i];
      Pdot[i] = -0.5*(kiminus1*(kiminus1-1) - ki*(ki-1))/theta;
    }
  }

  double density() {

    double logP = 0.0;

    sortNodeList(nodeList);
    ArrayList<Node> sorted = getNodeList();

    for (int i=0; i<sorted.size ()-1; i++) {
      int k = sorted.get(i).getLineages();
      double delta = sorted.get(i+1).getAge()-sorted.get(i).getAge();
      logP += 0.5*k*(k-1)*delta/theta;
    }

    return logP;
  }

  void doStep() {

    double[] Q = new double[nodeList.size()];
    double[] P = new double[nodeList.size()];
    double[] Qdot = new double[nodeList.size()];
    double[] Pdot = new double[nodeList.size()];

    getQandP(Q, P);
    getQdotAndPdot(Q, P, Qdot, Pdot);

    for (int i=0; i<nodeList.size (); i++) {
      Q[i] = Q[i] + dt*Qdot[i];
      P[i] = P[i] + dt*Pdot[i];
    }

    setQandP(Q, P);

    updateTopology();
  }

  // What to do when interval sizes become negative
  void updateTopology() {

    for (Node nodeA : nodeList) {
      if (nodeA.isRoot())
        continue;

      Node nodeB = nodeA.getParent();

      if (nodeB.getAge()<nodeA.getAge()) {
        double delta = nodeA.getAge()-nodeB.getAge();
        if (nodeA.isLeaf()) {
          nodeB.setAge(nodeA.getAge()+delta); // reflection
          nodeB.setMomentum(-nodeB.getMomentum());
        } else {
          int idx = int(random(2));
          Node nodeC = nodeA.getChildren().get(idx);

          if (!nodeB.isRoot()) {
            Node nodeD = nodeB.getParent();
            nodeD.getChildren().remove(nodeB);
            nodeD.addChild(nodeA);
          } else {
            nodeA.setParent(null);
            setRoot(nodeA);
          }

          nodeB.getChildren().remove(nodeA);
          nodeA.getChildren().remove(nodeC);
          nodeA.addChild(nodeB);
          nodeB.addChild(nodeC);
          
        }
      }
    }
  }
}

