#ifndef LIMITEDSIZEPRIORITYQUEUE_H
#define LIMITEDSIZEPRIORITYQUEUE_H

#include <set>
#include <stddef.h>

/**
 * Based on https://www.geeksforgeeks.org/double-ended-priority-queue/
 */
template <typename ElementType>
class LimitedSizePriorityQueue {
public:
    LimitedSizePriorityQueue(size_t maxSize);

    void push(ElementType element);

    ElementType popSmall();
    ElementType popLarge();

    ElementType peekSmall();
    ElementType peekLarge();

    size_t size();
    bool empty();

private:
    std::set<ElementType> data;
    size_t maxSize;
};

template <typename ElementType>
void LimitedSizePriorityQueue<ElementType>::push(ElementType element)
{
    data.insert(element);
    if (size() > maxSize) {
        popLarge();
    }
}

template <typename ElementType>
ElementType LimitedSizePriorityQueue<ElementType>::popSmall()
{
    ElementType element = *(data.begin());
    if (size() > 0) {
        data.erase(data.begin());
    }
    return element;
}

template <typename ElementType>
ElementType LimitedSizePriorityQueue<ElementType>::popLarge()
{
    ElementType element = *(data.rbegin());
    if (size() > 0) {
        auto end = data.end();
        end--;
        data.erase(end);
    }
    return element;
}

template <typename ElementType>
ElementType LimitedSizePriorityQueue<ElementType>::peekSmall()
{
    return *(data.begin());
}

template <typename ElementType>
ElementType LimitedSizePriorityQueue<ElementType>::peekLarge()
{
    return *(data.rbegin());
}

template <typename ElementType>
size_t LimitedSizePriorityQueue<ElementType>::size()
{
    return data.size();
}

template <typename ElementType>
bool LimitedSizePriorityQueue<ElementType>::empty()
{
    return size() == 0;
}

#endif // LIMITEDSIZEPRIORITYQUEUE_H
