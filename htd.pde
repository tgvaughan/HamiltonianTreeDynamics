int margin = 50;
int nodeRadius = 10;
double scaleY = 200;
Tree tree;

void setup() {
     size(800, 600);

//     frameRate(frame_rate);

     stroke(#000000);
     fill(#ff0000);
     
     Node leaf1 = new Node().setAge(0);
     Node leaf2 = new Node().setAge(1);
     Node leaf3 = new Node().setAge(2);
     
     Node internal = new Node().setAge(1.5);
     internal.addChild(leaf1).addChild(leaf2);
     Node root = new Node().setAge(2.5);
     root.addChild(internal).addChild(leaf3);
     
     tree = new Tree(root);
     
     ArrayList<Node> sorted = tree.getSortedNodeList();
     for (int i=0; i<tree.getNodeCount(); i++) {
       println(sorted.get(i).getAge() + " " + sorted.get(i).getLineages());
     }
     
     println(density(tree, 1.0));
}

void draw() {
  background(#ffffff);
  drawTree();
}

void drawTree() {      
  
      for (Node node : tree.getNodeList()) {
        if (node.isRoot())
          continue;
        
        int[] pos = node.getScreenPosition();
        int[] parentPos = node.getParent().getScreenPosition();
        line(pos[0], pos[1], parentPos[0], parentPos[1]);
      }
  
      for (Node node : tree.getNodeList()) {
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
      
      Node setAge(double age) {
        this.age = age;
        return this;
      }
      
      double getAge() {
        return age;
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
          for (int i=0; i<children.size(); i++) {
            res.addAll(children.get(i).getChildList());
          }
        }
        
        return res;
      }
      
      ArrayList<Node> getNodeList() {
        ArrayList<Node> res = new ArrayList<Node>();
        
        res.add(this);
        for (int i=0; i<children.size(); i++) {
          res.addAll(children.get(i).getNodeList());
        }
        
        return res;
      }
      
      double getLayoutPosition() {
        if (layoutPosition<0) {
          layoutPosition = 0.0;
          for (int i=0; i<children.size(); i++) {
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
          for (int i=0; i<children.size(); i++) {
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
      
      Tree(Node root) {
        this.root = root;
        childList = root.getChildList();
        nodeList = root.getNodeList();
        
        for (int i=0; i<childList.size(); i++) {
          childList.get(i).layoutPosition = i;
        }
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
      
      ArrayList<Node> getSortedNodeList() {
        boolean sorted;
        do {
          sorted = true;
          for (int i=0; i<nodeList.size()-1; i++) {
            if (nodeList.get(i+1).getAge()<nodeList.get(i).getAge()) {
              Node tmp = nodeList.get(i+1);
              nodeList.set(i+1, nodeList.get(i));
              nodeList.set(i, tmp);
              sorted = false;
            }
          }
        } while(!sorted);
        
        int k=0;
        for (Node node : nodeList) {
          if (node.isLeaf())
            k += 1;
          else
            k -= 1;
          node.setLineages(k);
        }
        
        return nodeList;
      }
      
}
ff
double density(Tree tree, double theta) {
  
  double logP = 0.0;
  
  ArrayList<Node> sorted = tree.getSortedNodeList();
  
  for (int i=0; i<sorted.size()-1; i++) {
    int k = sorted.get(i).getLineages();
    double delta = sorted.get(i+1).getAge()-sorted.get(i).getAge();
    logP += 0.5*k*(k-1)*delta/theta;
  }
  
  return logP;
}
