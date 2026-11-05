"""Compact terminal progress, with per-item details written to a local log."""
from __future__ import annotations

from datetime import datetime, timezone
import json
from pathlib import Path

from rich.console import Console
from rich.progress import (
    BarColumn, MofNCompleteColumn, Progress, ProgressColumn,
    SpinnerColumn, TaskProgressColumn, TextColumn, TimeRemainingColumn,
)
from rich.text import Text


class RemainingColumn(ProgressColumn):
    def render(self, task):
        remaining = max(0, int((task.total or 0) - task.completed))
        return Text(f"{remaining:,} left", style="dim")


class RunProgress:
    def __init__(self, console: Console):
        self.console = console
        self.interactive = console.is_interactive
        self.progress = Progress(
            SpinnerColumn(), TextColumn('{task.description}', markup=False),
            BarColumn(bar_width=None), TaskProgressColumn(), MofNCompleteColumn(),
            RemainingColumn(), TextColumn('ETA'), TimeRemainingColumn(),
            console=console, transient=True, disable=not self.interactive,
            refresh_per_second=4,
        )
        self.log_path: Path | None = None
        self._stream = None
        self._depth = 0

    def __enter__(self):
        if self._depth == 0:
            if self.log_path:
                self.log_path.parent.mkdir(parents=True, exist_ok=True)
                self._stream = self.log_path.open('a', encoding='utf-8')
            self.progress.start()
        self._depth += 1
        return self

    def __exit__(self, exc_type, exc, traceback):
        self._depth -= 1
        if self._depth == 0:
            self.progress.stop()
            if self._stream:
                self._stream.close()
                self._stream = None

    def log(self, message, **fields):
        if self._stream:
            entry = {'time':datetime.now(timezone.utc).isoformat(), 'message':message, **fields}
            self._stream.write(json.dumps(entry, ensure_ascii=False) + '\n')
            self._stream.flush()

    def track(self, items, label, *, completed=0):
        """Count attempts completed in this session plus previously saved work.

        A failed attempt counts as processed; its retry remains in SQLite.
        An interrupted item does not advance the bar. ETA applies to this phase.
        """
        total = completed + len(items)
        if not items:
            return
        with self:
            task_id = self.progress.add_task(label, total=total, completed=completed)
            if not self.interactive:
                self.console.print(f'{label}: {completed:,}/{total:,} ({len(items):,} remaining)', markup=False)
            finished = False
            try:
                for item in items:
                    yield item
                    self.progress.advance(task_id)
                finished = True
            finally:
                task = self.progress.tasks[0] if len(self.progress.tasks) == 1 else next(
                    t for t in self.progress.tasks if t.id == task_id)
                self.log('phase_finished' if finished else 'phase_interrupted',
                         phase=label, completed=int(task.completed), total=total)
                if not self.interactive:
                    state = 'complete' if finished else 'interrupted'
                    self.console.print(f'{label}: {int(task.completed):,}/{total:,} {state}', markup=False)
                self.progress.remove_task(task_id)
