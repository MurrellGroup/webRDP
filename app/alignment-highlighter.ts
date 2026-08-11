export type ParentAffinityKind = "unique" | "shared" | "novel" | "missing";

export interface ParentAffinity {
  kind: ParentAffinityKind;
  parentSlots: number[];
}

export function isCanonicalBase(base: string | undefined): boolean {
  return base === "A" || base === "C" || base === "G" || base === "T";
}

export function classifyParentAffinity(
  base: string | undefined,
  parentBases: Array<string | undefined>,
): ParentAffinity {
  if (!isCanonicalBase(base)) return { kind: "missing", parentSlots: [] };
  const canonicalParents = parentBases
    .map((parentBase, parentSlot) => ({ parentBase, parentSlot }))
    .filter(({ parentBase }) => isCanonicalBase(parentBase));
  if (!canonicalParents.length) return { kind: "missing", parentSlots: [] };
  const parentSlots = canonicalParents
    .filter(({ parentBase }) => parentBase === base)
    .map(({ parentSlot }) => parentSlot);
  if (parentSlots.length === 0) return { kind: "novel", parentSlots: [] };
  if (parentSlots.length === 1) return { kind: "unique", parentSlots };
  return { kind: "shared", parentSlots };
}

export function parentInformativeSites(
  sequences: string[],
  parentIndexes: number[],
  length: number,
): number[] {
  const indexes: number[] = [];
  for (let site = 0; site < length; site += 1) {
    const alleles = new Set<string>();
    for (const parentIndex of parentIndexes) {
      const base = sequences[parentIndex]?.[site];
      if (isCanonicalBase(base)) alleles.add(base as string);
    }
    if (alleles.size > 1) indexes.push(site);
  }
  return indexes;
}

export function affinityDescription(
  affinity: ParentAffinity,
  parentNames: string[],
): string {
  if (affinity.kind === "missing") return "gap, ambiguity, or no callable parent base";
  if (affinity.kind === "novel") return "matches none of the selected parents";
  const names = affinity.parentSlots.map((slot) => parentNames[slot] ?? `Parent ${slot + 1}`);
  return affinity.kind === "shared"
    ? `shared by ${names.join(", ")}`
    : `matches ${names[0]}`;
}
