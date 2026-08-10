import assert from "node:assert/strict";
import {
  DEFAULT_OPTIONS,
  demoEvent,
  eventLength,
  exportRecombinationFree,
  makeDemoAlignment,
  parseGenomeAnnotations,
  parseAlignment,
  parseProject,
  serializeProject,
} from "../app/rdp-core";

const event = {
  ...demoEvent(),
  start: 2_200,
  end: 100,
  wraps: true,
  breakpointModel: {
    method: "two-state-hmm" as const,
    informativeSites: 284,
    stateSwitches: 2,
    majorFit: 0.96,
    minorFit: 0.94,
  },
  decision: "accepted" as const,
  note: "reviewed circular positive control",
};
const serialized = serializeProject({
  alignment: {
    ...makeDemoAlignment(),
    features: parseGenomeAnnotations("genome\ttest\tCDS\t101\t900\t.\t+\t0\tID=cds1;Name=rep", "fixture.gff3", 2_400),
  },
  options: { ...DEFAULT_OPTIONS, candidateParents: 12 },
  events: [event],
  metrics: { elapsedMs: 12.5, comparisons: 56, engine: "test" },
  distance: [0, 0.1, 0.1, 0],
  auditLog: [{ id: "audit-1", timestamp: "2026-08-11T00:00:00.000Z", action: "Accepted event", summary: "Reviewed fixture.", eventId: event.id }],
});
const restored = parseProject(serialized);

assert.equal(restored.schema, "rdp-web/0.5");
assert.equal(restored.alignment.sequences[0].sequence, makeDemoAlignment().sequences[0].sequence);
assert.equal(restored.options.candidateParents, 12);
assert.equal(restored.events[0].decision, "accepted");
assert.equal(restored.events[0].note, "reviewed circular positive control");
assert.equal(restored.events[0].wraps, true);
assert.equal(restored.events[0].breakpointModel?.method, "two-state-hmm");
assert.equal(restored.events[0].history[0].action, event.history[0].action);
assert.equal(restored.events[0].evidenceStale, false);
assert.equal(restored.alignment.features?.[0].name, "rep");
assert.equal(eventLength(restored.events[0], restored.alignment.length), 300);
assert.equal(restored.events[0].evidence[0].calibration, event.evidence[0].calibration);
assert.deepEqual(restored.distance, [0, 0.1, 0.1, 0]);
assert.equal(restored.auditLog[0].eventId, event.id);

const masked = exportRecombinationFree(restored.alignment, restored.events, "mask")[0];
const maskedAlignment = parseAlignment(masked.content, masked.filename);
assert.equal(maskedAlignment.sequences[0].sequence.match(/N/g)?.length, 300);
assert.equal(maskedAlignment.sequences[1].sequence.match(/N/g)?.length ?? 0, 0);
