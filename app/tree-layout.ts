import type { NeighborJoiningNode } from "./rdp-core";

export interface TreeLayoutEdge {
  parentX: number;
  childX: number;
  childY: number;
}

export interface TreeLayoutJoint {
  x: number;
  y: number;
  minimumChildY: number;
  maximumChildY: number;
}

export interface TreeLayoutLabel {
  x: number;
  y: number;
  name: string;
}

export interface TreeLayout {
  width: number;
  height: number;
  path: string;
  edges: TreeLayoutEdge[];
  joints: TreeLayoutJoint[];
  labels: TreeLayoutLabel[];
  zeroLengthBranches: number;
}

function safeLength(value: number): number {
  return Number.isFinite(value) ? Math.max(0, value) : 0;
}

function leafCount(node: NeighborJoiningNode): number {
  return node.children?.length
    ? node.children.reduce((total, child) => total + leafCount(child), 0)
    : 1;
}

function maximumDistance(node: NeighborJoiningNode, parentDistance = 0): number {
  const distance = parentDistance + safeLength(node.length);
  return node.children?.length
    ? Math.max(distance, ...node.children.map((child) => maximumDistance(child, distance)))
    : distance;
}

function maximumDepth(node: NeighborJoiningNode, depth = 0): number {
  return node.children?.length
    ? Math.max(depth, ...node.children.map((child) => maximumDepth(child, depth + 1)))
    : depth;
}

export function layoutNeighborJoiningTree(
  root: NeighborJoiningNode,
  width = 720,
  labelWidth = 190,
): TreeLayout {
  const leaves = leafCount(root);
  const height = Math.max(150, leaves * 24 + 28);
  const left = 20;
  const drawableWidth = Math.max(80, width - left - labelWidth);
  const farthest = maximumDistance(root);
  const deepest = Math.max(1, maximumDepth(root));
  const useDistance = farthest > 1e-12;
  const scale = drawableWidth / (useDistance ? farthest : deepest);
  const edges: TreeLayoutEdge[] = [];
  const joints: TreeLayoutJoint[] = [];
  const labels: TreeLayoutLabel[] = [];
  let nextLeaf = 0;
  let zeroLengthBranches = 0;

  function place(
    node: NeighborJoiningNode,
    parentDistance: number,
    depth: number,
  ): { x: number; y: number } {
    const length = safeLength(node.length);
    if (depth > 0 && length <= 1e-12) zeroLengthBranches += 1;
    const distance = parentDistance + length;
    const x = left + (useDistance ? distance : depth) * scale;
    if (!node.children?.length) {
      const y = 16 + (nextLeaf + 0.5) * ((height - 32) / Math.max(1, leaves));
      nextLeaf += 1;
      labels.push({ x: Math.min(width - labelWidth + 8, x + 7), y, name: node.name ?? "unnamed" });
      return { x, y };
    }
    const children = node.children.map((child) => place(child, distance, depth + 1));
    const minimumChildY = Math.min(...children.map((child) => child.y));
    const maximumChildY = Math.max(...children.map((child) => child.y));
    const y = children.reduce((sum, child) => sum + child.y, 0) / children.length;
    joints.push({ x, y, minimumChildY, maximumChildY });
    children.forEach((child) => edges.push({ parentX: x, childX: child.x, childY: child.y }));
    return { x, y };
  }

  const rootPosition = place(root, 0, 0);
  const pathParts = joints.map((joint) => `M${joint.x.toFixed(2)},${joint.minimumChildY.toFixed(2)}V${joint.maximumChildY.toFixed(2)}`);
  edges.forEach((edge) => pathParts.push(`M${edge.parentX.toFixed(2)},${edge.childY.toFixed(2)}H${edge.childX.toFixed(2)}`));
  pathParts.push(`M${left.toFixed(2)},${rootPosition.y.toFixed(2)}H${rootPosition.x.toFixed(2)}`);

  return {
    width,
    height,
    path: pathParts.join(" "),
    edges,
    joints,
    labels,
    zeroLengthBranches,
  };
}
