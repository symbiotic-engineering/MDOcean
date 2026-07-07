find "$1" -type f -name '*.json' -print0 |
xargs -0 jq -r '
    .. | objects
    | select(has("requests"))
    | .requests[]
    | .modelId?
' 2>/dev/null |
sort -u