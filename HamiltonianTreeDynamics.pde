int margin = 50;
int nodeRadius = 10;
double scaleY = 150;
int frame_rate = 25;

Tree tree;
double theta = 1.0;

double dt = 1e-2;

void setup() {
  size(800, 600);
  frameRate(frame_rate);

  stroke(#000000);
  fill(#ff0000);

  Node leaf1 = new Node().setAge(0);
  Node leaf2 = new Node().setAge(1);
  Node leaf3 = new Node().setAge(2);

  Node internal = new Node().setAge(3);
  internal.addChild(leaf1).addChild(leaf2);
  Node root = new Node().setAge(3.2);
  root.addChild(internal).addChild(leaf3);

  tree = new Tree(root);
}

void draw() {
  tree.doStep();
  background(#ffffff);
  drawTree();
}

//void keyPressed() {
//  tree.doStep();
//}

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

  void sortNodeList() {
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
    sortNodeList();
    for (int i=0; i<nodeList.size(); i++) {
      Q[i] = nodeList.get(i).getAge();
      P[i] = nodeList.get(i).getMomentum();
    }
  }
  
  void setQandP(double[] Q, double[] P) {
    for (int i=0; i<nodeList.size(); i++) {
      nodeList.get(i).setAge(Q[i]);
      nodeList.get(i).setMomentum(P[i]);
    }
  }
  
  void getQdotAndPdot(double[] Q, double[] P, double[] Qdot, double[] Pdot) {
    for (int i=0; i<nodeList.size(); i++) {
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

    sortNodeList();
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
    getQdotAndPdot(Q,P,Qdot,Pdot);
    
    for (int i=0; i<nodeList.size(); i++) {
      print(Qdot[i] + " ");
      Q[i] = Q[i] + dt*Qdot[i];
      P[i] = P[i] + dt*Pdot[i];
    }
    println();
    
    setQandP(Q, P);
    
    updateTopology();
  }
  
  void updateTopology() {
    
    for (Node node : nodeList) {
      if (node.isRoot())
        continue;
      
      Node parent = node.getParent();
      
      if (parent.getAge()<node.getAge()) {
        double delta = node.getAge()-parent.getAge();
        if (node.isLeaf()) {
          parent.setAge(node.getAge()+delta); // reflection
          parent.setMomentum(-parent.getMomentum());
        } else {
          int idx = int(random(2));
          Node child = node.getChildren().get(idx);
          
          if (!parent.isRoot()) {
            Node grandParent = parent.getParent();
            grandParent.getChildren().remove(parent);
            grandParent.addChild(node);
          }
          
          parent.getChildren().remove(node);
          node.getChildren().remove(child);
          node.addChild(parent);
          parent.addChild(child);
        }
      }
    }
    
  }

}

