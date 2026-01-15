---
layout: page
title: Schedule
nav_order: 2
description: Schedule for the entire semester.
---

## Schedule

{% for schedule in site.schedules %}
{% if schedule.module_title %}
<h2 class="fs-4 mt-0 module-title">{{ schedule.module_title }}</h2>
{% endif %}
<h3 class="fs-5 mt-2">{{ schedule.title }}</h3>
<div class="module">
  {{ schedule.content }}
</div>
{% endfor %}
